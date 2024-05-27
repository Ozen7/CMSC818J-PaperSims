import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from collections import deque
import math
import time


class Scheduler:
    # The scheduler is where the Huffman Tree is located. When the MergeTree finishes its merging, it calculates the next set of partial matrices to merge
    # it then goes to the PartialMatrixFetcher to pull any partial matrices that are small enough into the mergeTree, and to the MatrixAFetcher to pull 
    # a set of compressed columns (each of which becomes a partial matrix to be sent to the MergeTree). 
    def __init__(self, MAF: "MatrixAFetcher", MT:"MergeTree", DLB: "DistanceListBuilder", MBP: "MatrixBPrefetcher", MA: "MultiplierArray", numMergers) -> None:
        self.MatrixAFetcher = MAF
        self.MergeTree = MT
        self.DistanceListBuilder = DLB
        self.MatrixBPrefetcher = MBP
        self.MultiplierArray = MA
        self.MatrixARows = MAF.rowLengths.copy()
        self.StoredMatrices = []
        self.maxSchedule = ((len(self.MatrixARows)-2) % (numMergers-1)) + 2
        self.numMergers = numMergers
        self.merging = False
        self.rounds = 0

    def addMatrix(self,num,size):
        self.StoredMatrices.append((num,size))
        self.StoredMatrices.sort(key= lambda x: x[1],reverse=True)
        
        
    # Called by the MergeTree when all of its input FIFOs are empty, and the highest level fifo is also empty
    def schedule(self):
        numACols = 0
        matrices = []
        c = 0
        PMlength = 0
        while c < self.maxSchedule and self.MatrixARows and self.StoredMatrices:
            if self.MatrixARows[-1] < self.StoredMatrices[-1][1]:
                numACols += 1
                PMlength += len(self.matrixARows[-1])
                self.MatrixARows.pop()
            else:
                matrices.append(self.StoredMatrices[-1][0])
                self.StoredMatrices.pop()
            c += 1
        
        while c < self.maxSchedule and self.MatrixARows:
            numACols += 1
            self.MatrixARows.pop()
            c += 1
        while c < self.maxSchedule and self.StoredMatrices:
            matrices.append(self.StoredMatrices[-1][0])
            self.StoredMatrices.pop()
            c += 1
        
        self.rounds += 1
        self.maxSchedule = self.numMergers
        self.MatrixAFetcher.inputInstructions(numACols)
        self.MergeTree.inputPartials(numACols,matrices,PMlength)
        self.merging = False
        
        
    # Called by the MultiplierArray when it has no input
    def checkEOF(self):
        # print(self.merging, self.MatrixAFetcher.inputFlag, self.DistanceListBuilder.lookaheadFIFO, self.MatrixBPrefetcher.inputFlag)
        # if not merging, and no input remaining:
        if self.MatrixAFetcher.endFlag and not self.merging and not self.DistanceListBuilder.lookaheadFIFO and not self.MatrixBPrefetcher.inputFlag:
            self.MatrixBPrefetcher.endFlag = True
            self.MergeTree.endFlag = True
            self.MultiplierArray.endFlag = True
        if not self.merging and not self.MatrixAFetcher.inputFlag and not self.DistanceListBuilder.lookaheadFIFO and not self.MatrixBPrefetcher.inputFlag:
            self.merging = True
            # indicates to the merge Tree that it should begin merging
            
            self.MergeTree.startMerging()
        
        
        

class Memory:
    def __init__(self, numChannels, peakBandwidthPerChannel) -> None:
        self.numChannels = numChannels
        self.channelQueues = [deque() for _ in range(numChannels)]
        self.channelCurrent = [(None,0,0) for _ in range(numChannels)]
        self.BPC = peakBandwidthPerChannel
        self.TotalMemoryPulled = 0
        self.NumCyclesInUse = 0 #use this for memory efficiency when in use.
    
    def requestData(self,loader, channel, amount):
        # numCycles includes the row and column latencies
        self.channelQueues[channel].append([loader,amount]) # (loader, number of bytes, memory latency)
    
    def cycle(self):
        inuse = False
        for i,channel in enumerate(self.channelCurrent):
            if channel[1] > 0:
                inuse = True
                sent = min(channel[1],self.BPC)
                self.TotalMemoryPulled += sent
                channel[0].freeBytes += sent # "send" the bytes over
                channel[1] -= sent
            elif channel[1] <= 0 and self.channelQueues[i]:
                self.channelCurrent[i] = self.channelQueues[i].popleft()
        if inuse:
            self.NumCyclesInUse += 1
            
        
        

# A huge issue with the matrix A fetcher in this accelerator is that how it actually functions is NOT mentioned at all anywhere in the paper!
# Due to this, I have decided that the fetcher can check one row every cycle to see if it contains a matrix within the condensed columns. If it does, 
# the up to 64 condensed columns in that row are put into the lookahead fifo in one cycle (as stated in the Architectural Setup table) and the row is incremented

# This way of computing might seem incredibly inefficient, especially since we could be working with matrices in the billions of parameters, but we should note that
# this can actually be done while waiting for merging to occur. Since the size of the new matrix can be calculated beforehand (rough estimate), the huffman scheduler can do its work
# as the merging is processed. Ideally, we fill up the queue while merging occurs, and stop until merging ends before continnuing. Then, as the data continues to flow through the accelerator,
# we can continue going through rows and finding matches. Also, the waiting time (when the FIFO isnt full) will be slightly covered up by the time it takes 
# the B matrix loader to load a row from memory.

# I might also add functionality to look at the next n rows (like 64-128) and skip all the ones that have a gap less than the minimum size of the condensed columns,because 
# they can't possibly have enough values to be in the condensed columns we want. This mainly protects against an immense waste in a hypersparse matrix.

class MatrixAFetcher:
    def __init__(self,A: csr_matrix) -> None:
        self.A = A
        # Note: for this accelerator, the endFlag will ONLY be set by the software scheduler
        self.endFlag = False
        self.inputFlag = False
        self.distanceListBuilder = None
        self.colLengths = []
        self.rowLengths = []
        self.memory = None

        
        # processing to make sure that we have a list of the number of values in each compressed column - important for the huffman tree scheduler
        # while doing this, we also create a list of nonzero row lengths to be used when iterating through them to find A values to send.
        for x in range(len(A.indptr)-1):
            # If the paper had mentioned that this sort of processing was done beforehand, it would be fine
            # my only issue is that they didn't mention it, even if it was necessary
            l = A.indptr[x+1] - A.indptr[x]
            self.rowLengths.append(l)
            if len(self.colLengths) < l:
                self.colLengths += ([0] * (l - len(self.colLengths)))
            for y in range(l):
                self.colLengths[y] += 1
        # The currentCondensedColumn should be equal to the furthest condensed column at first
        self.currentCondensedColumn = len(self.colLengths)
        self.currentNumColumns = 0
        self.currentRow = 0
        
        self.memoryUse = 0
        self.memoryQueue = deque([])
        self.freeBytes = 0

    def setDLB(self,DLB):
        self.distanceListBuilder = DLB
    
    def setMemory(self,Memory):
        self.memory = Memory
    
    def inputInstructions(self,numColumns):
        if self.currentCondensedColumn == 0:
            self.currentNumColumns = 0
            self.endFlag = True
            return
        self.inputFlag = True
        self.currentRow = 0
        self.currentCondensedColumn -= numColumns
        self.currentNumColumns = numColumns
        if self.currentCondensedColumn < 0:
            self.currentNumColumns += self.currentCondensedColumn
            self.currentCondensedColumn = 0
            
        # The paper simply states that the matrix A fetcher "calculates addresses of data in the selected columns", and not how it is 
        # done.

        # Make use of current condensed column number, and set currentCondensedColumn (the minimum number of columns needed in a row to have values we want for this round)
        # currentNumColumns allows us to specify the range of columns we want to pull from memory (condensed columns [currendCondensedColumns:currentCondensedColumns + currentNumColumns])
        # This number is specified by the huffman tree builder, making use of the number of values in each column within a minheap
    
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def cycle(self):
        # Note: MatrixAFetcher does not wait for values to be fetched, and keeps making requests
        while self.freeBytes >= 8:
            self.distanceListBuilder.input(self.memoryQueue.popleft()) # data, row, column, condensed column number are sent to DLB buffer. 
            self.freeBytes -= 8

        if self.endFlag or not self.inputFlag:
            return
        
        # Every cycle, if we arent at EOF:
        # 1) Look through a row, check if it has have columns within the condensed column.
        # 2) If it does, load them into the FIFO and increment
        # 3) If it doesn't, increment.
        if self.currentRow >= len(self.A.indptr) -1 and not self.memoryQueue:
            # We have gone through all the rows, and all valid A values are sent. Send none so that the Multiplier can signal that loading is finished to the MergeTree
            self.inputFlag = False
            return
        elif self.currentRow >= len(self.A.indptr) -1:
            return
        
        if self.rowLengths[self.currentRow] > self.currentCondensedColumn:
            # i = the exact point in the data and column arrays where the condensed column set starts
            i = self.A.indptr[self.currentRow] + self.currentCondensedColumn
            if self.rowLengths[self.currentRow] - self.currentCondensedColumn < self.currentNumColumns: 
                end = self.A.indptr[self.currentRow+1]
            else:
                end = self.currentNumColumns + self.A.indptr[self.currentRow] + self.currentCondensedColumn
                
            # because channels are assigned cyclically, and each of 16 channels is 64 values wide, after 16 * 64 values, we'll be in the same spot within the same channel
            # C finds us which channel we start at, and o finds us the offset from the beginning of the channel
            O = (i%64)
            C = (i//64)%16
            for x in range(i,end):
                self.memoryUse += 4 * 2 # each int is 32 bits (or 4 bytes), and there are 2 of them to fetch (data & col). 8 bytes are fetched
                # from a random channel. To simplify things, we assume that every 8 bytes fetched correspond
                # to a set of data to be sent to the Distance List Builder, and send a value every 8 bytes that are fetched from memory.

                # we request chunks of 8 bytes
                self.memory.requestData(self,C,8)
                # increment offset. If we move to a new channel, increment the channel number mod 8
                O += 8
                O %= 64
                if O < 8:
                    C += 1
                    C %= 16
                    
                self.memoryQueue.append((self.A.data[x],self.currentRow,self.A.indices[x],x-i)) # data has to be sent in order
        self.currentRow += 1
    
    
class DistanceListBuilder:
    # Output from MatrixAFetcher comes into here.
    # Look at the buffer, and pass the required rows as well as the distance of that row to the B prefetcher once it comes into the buffer
    # Each cycle, wait for memory to be fetched by B prefetcher, and then send the data to the multiplier array (in the same cycle as the B laoder)
    # Basically just have a flag that causes it to send and pop the buffer, that is set by the B prefetcher
    def __init__(self) -> None:
        self.lookaheadFIFO = deque([])
        self.MBP = None
        self.Multiplier = None
    
    def setMBP(self, MBP: "MatrixBPrefetcher"):
        self.MBP = MBP
    
    def setMultiplier(self, M: "MultiplierArray"):
        self.Multiplier = M
        
    def input(self,Data):
        self.lookaheadFIFO.append(Data)
        # send the required row as well as the current length of the buffer (how many A values are sent to the multiplier until this row is needed)
        # Send the column number of the A value, which corresponds to the row number of the B matrix that will be multiplied by it, as well as what its position is to know when it will be used
        self.MBP.input(Data[2], len(self.lookaheadFIFO)) 
        
    def sendInput(self):
        self.Multiplier.loadAvalue(self.lookaheadFIFO.popleft())
   
   
class MatrixBCache:
    # Instead of keeping track of future accesses, we only keep track of past accesses and decide which value to remove based on age
    # the older a row is, the less valuable.
    # Make a version where we do prefetch the next value, but don't look ahead, and make a version where we don't prefetch or look ahead.
    def __init__(self, isPrefetch, B:csr_matrix) -> None:
        self.endFlag = False
        self.isPrefetch = isPrefetch
        self.B = B
        
        self.timeStep = 0
        self.inputFlag = False
        self.prefetchBufferHash = {-1: 1024} # initially, every single set in the row is empty (assigned to -1), there are 1024 lines
        self.rowAgeHash = {-1: 0}
        self.multiplierArray = None
        self.memory = None
        
        self.inputBuffer = deque([])
        self.oldestRowInCache = -1

        
        self.currentRowPointer = 0 # Used for row that is being sent to multiplier 
        self.currentRowSendingFinished = True # indicates that the current row has been fully sent. Starts as true, but is set to false in the first cycle
        self.currentRowNumber = -1 # current row number. Used for accessing rowAgeHash and prefetchBufferHash
        self.currentRowData = [] # The current row's Data
        self.currentRowIndices = [] # The current row's column indices
        
        
        if self.isPrefetch:
            self.nextRowPointer = 0 # Used for row that is being loaded into buffer
            self.nextRowLoadingFinished = True # same as above, indicates that next row has been fully loaded.
            self.nextRowNumber = -1 # next row number. Used for accessing rowAgeHash and prefetchBufferHash
            self.nextRowData = [] # The next row's Data
            self.nextRowIndices = [] # The next row's column indices
        else:
            self.currentRowLoadingFinished = True
        
        self.memoryWastedCycles = 0
        self.memoryAccessBytes = 0
        self.freeBytes = 0
        self.rowsToFetch = 0
        
        self.wastedCycles = 0
        pass  
    
    
    def takeTimeStep(self):
        self.inputFlag = False
        self.timeStep += 1
        
        # All data from next row is now the current row, and we begin sending it to the multiplier each cycle
        if self.isPrefetch:
            self.currentRowPointer = 0
            self.currentRowData = self.nextRowData
            self.currentRowIndices = self.nextRowIndices
            self.currentRowNumber = self.nextRowNumber
            self.currentRowSendingFinished = self.currentRowNumber == -1 # if our row number is -1, we are "done" sending data. O.W, we are sending data

            # Indicate to DLB that a new A value should be sent, since we're done sending everything for our row now.
            if self.currentRowNumber != -1:
                self.DistanceListBuilder.sendInput()
            
            # If there is nothing left in the distance fifo, nothing is loaded for this timestep,
            # even if something is sent in the middle of sending to multiplier it won't be loaded until the next timestep
            self.nextRowNumber = -1
            self.nextRowData = []
            self.nextRowIndices = []
            self.nextRowPointer = 0
            if self.inputBuffer:
                self.nextRowNumber = self.inputBuffer.popleft()
                self.nextRowLoadingFinished = False
                self.nextRowIndices = list(self.B.indices[self.B.indptr[self.nextRowNumber]:self.B.indptr[self.nextRowNumber+1]])
                self.nextRowData = list(self.B.data[self.B.indptr[self.nextRowNumber]:self.B.indptr[self.nextRowNumber+1]])
                if self.nextRowNumber in self.prefetchBufferHash:
                    self.nextRowPointer = self.prefetchBufferHash[self.nextRowNumber]*48 #default is 0
                    
            # If something is being loaded or sent:
            if not (self.nextRowLoadingFinished and self.currentRowSendingFinished):
                self.inputFlag = True
        else:


            self.currentRowNumber = -1
            self.currentRowData = []
            self.currentRowIndices = []
            self.currentRowPointer = 0
            self.currentRowLoadingPointer = 0
            if self.inputBuffer:
                self.inputFlag = True
                self.currentRowNumber = self.inputBuffer.popleft()
                self.currentRowLoadingFinished = False
                self.currentRowIndices = list(self.B.indices[self.B.indptr[self.currentRowNumber]:self.B.indptr[self.currentRowNumber+1]])
                self.currentRowData = list(self.B.data[self.B.indptr[self.currentRowNumber]:self.B.indptr[self.currentRowNumber+1]])
                if self.currentRowNumber in self.prefetchBufferHash:
                    self.currrentRowLoadingPointer = self.prefetchBufferHash[self.currentRowNumber]*48 #default is 0
            
            if self.currentRowNumber != -1:
                self.DistanceListBuilder.sendInput()
            self.currentRowSendingFinished = self.currentRowNumber == -1 

    def input(self, col, length):
        # add col to our input buffer
        self.inputBuffer.append(col)
        self.inputFlag = True
        
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    
    def setMultiplier(self,M:"MultiplierArray"):
        self.multiplierArray = M
    
    def setDistanceListBuilder(self,DLB:"DistanceListBuilder"):
        self.DistanceListBuilder = DLB
    
    def setMemory(self, M: "Memory"):
        self.memory = M
        
        
    def memoryAccess(self, numRows, rowNum, currentLocation):
        if numRows == 0:
            return
        
        # Each row is 48 elements, and each element is 12 bytes
        numBytes = numRows * 48 * 12 
        
        # C finds us which channel we start at, and o finds us the offset from the beginning of the channel
        # rowNum is the row we're pulling, and currentLocation is how far in we are at this point.
        locationInList = self.B.indptr[rowNum] + currentLocation
        O = (locationInList%64)
        C = (locationInList//64)%8
        # we request chunks of 8 bytes, that's the peak memory bandwidth anyway so we don't lose anything
        # might be tempting to do chunks of 12 bytes, but that runs the risk of splitting across channels 
        # and also makes the data end up streaming slower (taking 2 cycles instead of 1)
        for x in range(0,(numBytes//8) + 1):
            ms = min(8,numBytes)
            self.memory.requestData(self,C,ms)
            self.memoryAccessBytes += ms
            numBytes -= ms
            # increment offset. If we move to a new channel, increment the channel number mod 8
            O += 8
            O %= 64
            if O < 8:
                C += 1
                C %= 8
        self.rowsToFetch = numRows # when this hits 0, the B prefetcher is done loading all data
    
    
    def cycle(self):
        if self.endFlag or not self.inputFlag:
            self.wastedCycles += 1
            return
        
        # if we are doing prefetching, the code is basically identical with distance replaced with age
        if self.isPrefetch:
            if self.currentRowSendingFinished and self.nextRowLoadingFinished:
                self.takeTimeStep()    
            elif self.currentRowSendingFinished and not self.nextRowLoadingFinished:
                self.memoryWastedCycles += 1
            

            # This is the part where we load a set of values into the multiplier array. Values are loaded 16 at a time
            # if our row pointer has exceeded the current length of the data, we don't do anything
            if not self.currentRowSendingFinished:
                end = self.currentRowPointer + 16 if self.currentRowPointer + 16 < len(self.currentRowData) else len(self.currentRowData)
                self.multiplierArray.loadBrow(self.currentRowData[self.currentRowPointer:end], self.currentRowIndices[self.currentRowPointer:end])
                self.currentRowPointer += 16
                # if we surpass the length of the current row after sending this data, we pop the row's row distance value and set a boolean
                # representing that the loading is finished and we can move on to the next row if the next set is in the buffer.
                if self.currentRowPointer >= len(self.currentRowData):
                    self.rowAgeHash[self.currentRowNumber] = self.timeStep
                    self.currentRowSendingFinished = True
        
            # This is the part where we load a row into buffer
            # MemoryFlag is false when prefetcher is waiting for memory, and true when it isnt
            if not self.nextRowLoadingFinished and self.rowsToFetch <= 0:                
                if self.nextRowPointer >= len(self.nextRowData):
                    self.nextRowLoadingFinished = True
                    return
                # Checks if the furthest row still has values to replace
                if self.prefetchBufferHash[self.oldestRowInCache] != 0:                                                
                    
                    # Assumption made here: one element = one value/coordinate pair from the B matrix. 48 elements in each line.
                    amountDataNeeded = min(self.prefetchBufferHash[self.oldestRowInCache] * 48, len(self.nextRowData) - self.nextRowPointer)
                    
                    self.nextRowPointer += amountDataNeeded # increment the row pointer by the amount of elements loaded into memory
                    # Decrement the number of lines that are "eaten up" that were spilt
                    amountDataNeeded = math.ceil(amountDataNeeded/48) # this is now in terms of the number of rows needed to be fetched
                    self.prefetchBufferHash[self.oldestRowInCache] -= amountDataNeeded
                    # Increment the number of lines that "replace" the previous row that was there
                    if self.nextRowNumber in self.prefetchBufferHash:
                        self.prefetchBufferHash[self.nextRowNumber] += amountDataNeeded
                    else:
                        self.prefetchBufferHash[self.nextRowNumber] = amountDataNeeded
                    
                    self.memoryAccess(amountDataNeeded, self.nextRowNumber, self.nextRowPointer) # pull indices and data into the buffer, this is counted in terms of prefetch buffer lines
                else:
                    # this code runs when the fetchers need a new place to fill in the buffer, it selects a oldest row, and deletes the current
                    # oldest row in the prefetchBufferHash
                    self.prefetchBufferHash.pop(self.oldestRowInCache)
                    
                    # we find a minimum time step last used, not a maximum next time of use!
                    minimum = 1000000000000000
                    for key in self.prefetchBufferHash:
                        v = self.rowAgeHash[key]
                        if v[0] < self.rowAgeHash[minimum]:
                            minimum = key
                    
                    self.oldestRowInCache = minimum
            elif self.rowsToFetch > 0:
                # reduce the number of rows to fetch by one for every 48 bytes that get fetched
                # since we don't use the bytes until the entire row is fetched, the order doesn't really matter anyway
                while self.freeBytes >= 48 * 12:
                    self.rowsToFetch -= 1
                    self.freeBytes -= 48 * 12
                return
        else:
            # if we are not doing prefetching, the code needs to be revamped
            if self.currentRowSendingFinished and self.currentRowLoadingFinished:
                self.takeTimeStep()    
            elif not self.currentRowLoadingFinished:
                self.memoryWastedCycles += 1
            
            # we only start sending once the loading is finished!
            if not self.currentRowSendingFinished and self.currentRowLoadingFinished:
                end = self.currentRowPointer + 16 if self.currentRowPointer + 16 < len(self.currentRowData) else len(self.currentRowData)
                self.multiplierArray.loadBrow(self.currentRowData[self.currentRowPointer:end], self.currentRowIndices[self.currentRowPointer:end])
                self.currentRowPointer += 16
                if self.currentRowPointer >= len(self.currentRowData):
                    self.rowAgeHash[self.currentRowNumber] = self.timeStep
                    self.currentRowSendingFinished = True
            # if loading isn't finished and we still have work to do:
            elif not self.currentRowLoadingFinished and self.rowsToFetch <= 0:                
                if self.currentRowLoadingPointer >= len(self.currentRowData):
                    self.currentRowLoadingFinished = True
                    return
                # Checks if the furthest row still has values to replace
                if self.prefetchBufferHash[self.oldestRowInCache] != 0:                                                
                    
                    # Assumption made here: one element = one value/coordinate pair from the B matrix. 48 elements in each line.
                    amountDataNeeded = min(self.prefetchBufferHash[self.oldestRowInCache] * 48, len(self.currentRowData) - self.currentRowLoadingPointer)
                    
                    self.currentRowLoadingPointer += amountDataNeeded # increment the row pointer by the amount of elements loaded into memory
                    # Decrement the number of lines that are "eaten up" that were spilt
                    amountDataNeeded = math.ceil(amountDataNeeded/48) # this is now in terms of the number of rows needed to be fetched
                    self.prefetchBufferHash[self.oldestRowInCache] -= amountDataNeeded
                    # Increment the number of lines that "replace" the previous row that was there
                    if self.currentRowNumber in self.prefetchBufferHash:
                        self.prefetchBufferHash[self.currentRowNumber] += amountDataNeeded
                    else:
                        self.prefetchBufferHash[self.currentRowNumber] = amountDataNeeded
                    
                    self.memoryAccess(amountDataNeeded, self.currentRowNumber, self.currentRowLoadingPointer) # pull indices and data into the buffer, this is counted in terms of prefetch buffer lines
                else:
                    # this code runs when the fetchers need a new place to fill in the buffer, it selects a oldest row, and deletes the current
                    # oldest row in the prefetchBufferHash
                    self.prefetchBufferHash.pop(self.oldestRowInCache)
                    
                    # we find a minimum time step last used, not a maximum next time of use!
                    minimum = 1000000000000000
                    for key in self.prefetchBufferHash:
                        v = self.rowAgeHash[key]
                        if v[0] < self.rowAgeHash[minimum]:
                            minimum = key
                    
                    self.oldestRowInCache = minimum
            # if loading is in progress:
            elif self.rowsToFetch > 0:
                # reduce the number of rows to fetch by one for every 48 bytes that get fetched
                # since we don't use the bytes until the entire row is fetched, the order doesn't really matter anyway
                while self.freeBytes >= 48 * 12:
                    self.rowsToFetch -= 1
                    self.freeBytes -= 48 * 12
                return
            
                
             
            
                               
       

class MatrixBPrefetcher:
    
    """
    PrefetchBufferHash contains the parts of each row that is contained in the buffer.
    Basically, the 16 fetchers only serve to load a max of 16 parts of a row at a time
    We only move one row ahead, meaning that while one row is being fed into the multiplier, we're loading the next one in
    If we don't get it done in time or we don't have enough room (unlikely), we quickly load the rest of the row and then start sending it to 
    The multiplier as we start loading yer another one
    
    get it? The critical parts are the following:
    1) fifo to keep track of which row is next
    2) timeStep: Increments every time a row is finished sending to the multiplier. 
    Basically just a way of properly being able to compare times, since not every value enters the fifo at the same time
    3) row distance hash: keeps track of the set of row distances that a specific row has in the FIFO
    4) prefetchBufferHash: Basically just keeps track of which parts of the buffer belong to which row. Useful when loading in data and replacing
    
    steps:
    1) Select the next row to be pulled from memory
    2) We find the row with the highest row distance, we replace it with our new row using the prefetchBufferHash
    3) Repeat step 2 until the entire row we'll need next is in the buffer
    4) Wait for the entire current row to be sent to multiplier, and then:
        1) pop its distance from RowDistanceHash 
        2) increment timeStep
        3) (when the loading of the next row is done) Pop the next row from the FIFO, this is the next row to be pulled from memory.
    7) return to step 1
    """
    
    
    def __init__(self, B:csr_matrix) -> None:
        self.B = B
        self.timeStep = 0
        self.endFlag = False
        self.inputFlag = False
        self.distanceFIFO = deque([])
        self.prefetchBufferHash = {-1: 1024} # initially, every single set in the row is empty (assigned to -1), there are 1024 lines
        self.rowDistanceHash = {-1: deque([0])} # a queue of row distances in order assigned to this row
        self.multiplierArray = None
        self.memory = None
        
        
        self.currentRowPointer = 0 # Used for row that is being sent to multiplier 
        self.nextRowPointer = 0 # Used for row that is being loaded into buffer
    
        self.currentRowSendingFinished = True # indicates that the current row has been fully sent. Starts as true, but is set to false in the first cycle
        self.currentRowNumber = -1 # current row number. Used for accessing rowDistanceHash and prefetchBufferHash
        self.currentRowData = [] # The current row's Data
        self.currentRowIndices = [] # The current row's column indices
        
        self.nextRowLoadingFinished = True # same as above, indicates that current row has been fully loaded.
        self.nextRowNumber = -1 # next row number. Used for accessing rowDistanceHash and prefetchBufferHash
        self.furthestRowInPrefetcher = -1 # the row number that the next row will take from to fill up. 
        self.nextRowData = [] # The next row's Data
        self.nextRowIndices = [] # The next row's column indices
        
        self.memoryWastedCycles = 0
        self.memoryAccessBytes = 0
        self.freeBytes = 0
        self.rowsToFetch = 0
        
        self.wastedCycles = 0
        

    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def takeTimeStep(self):
        self.inputFlag = False
        self.timeStep += 1
        
        # All data from next row is now the current row, and we begin sending it to the multiplier each cycle
        self.currentRowPointer = 0
        self.currentRowData = self.nextRowData
        self.currentRowIndices = self.nextRowIndices
        self.currentRowNumber = self.nextRowNumber
        self.currentRowSendingFinished = self.currentRowNumber == -1 # if our row number is -1, we are "done" sending data. O.W, we are sending data

        # Indicate to DLB that a new A value should be sent, since we're done sending everything for our row now.
        if self.currentRowNumber != -1:
            self.DistanceListBuilder.sendInput()
        
        # If there is nothing left in the distance fifo, nothing is loaded for this timestep,
        # even if something is sent in the middle of sending to multiplier it won't be loaded until the next timestep
        self.nextRowNumber = -1
        self.nextRowData = []
        self.nextRowIndices = []
        self.nextRowPointer = 0
        if self.distanceFIFO:
            self.nextRowNumber = self.distanceFIFO.popleft()
            self.nextRowLoadingFinished = False
            self.nextRowIndices = list(self.B.indices[self.B.indptr[self.nextRowNumber]:self.B.indptr[self.nextRowNumber+1]])
            self.nextRowData = list(self.B.data[self.B.indptr[self.nextRowNumber]:self.B.indptr[self.nextRowNumber+1]])
            if self.nextRowNumber in self.prefetchBufferHash:
                self.nextRowPointer = self.prefetchBufferHash[self.nextRowNumber]*48 #default is 0

            
        # If something is being loaded or sent:
        if not (self.nextRowLoadingFinished and self.currentRowSendingFinished):
            self.inputFlag = True
    
    def input(self, col, length):
        # add col to our distance fifo
        self.distanceFIFO.append(col)
        
        # add the new distance to our row distances
        if col in self.rowDistanceHash:
            self.rowDistanceHash[col].append(length)
        else:
            self.rowDistanceHash[col] = deque([length])
        
        # set inputFlag to true
        self.inputFlag = True
    
    def setMultiplier(self,M:"MultiplierArray"):
        self.multiplierArray = M
    
    def setDistanceListBuilder(self,DLB:"DistanceListBuilder"):
        self.DistanceListBuilder = DLB
    
    def setMemory(self, M: "Memory"):
        self.memory = M
        
        
    def memoryAccess(self, numRows, rowNum, currentLocation):
        if numRows == 0:
            return
        
        # Each row is 48 elements, and each element is 12 bytes
        numBytes = numRows * 48 * 12 
        
        # C finds us which channel we start at, and o finds us the offset from the beginning of the channel
        # rowNum is the row we're pulling, and currentLocation is how far in we are at this point.
        locationInList = self.B.indptr[rowNum] + currentLocation
        O = (locationInList%64)
        C = (locationInList//64)%8
        # we request chunks of 8 bytes, that's the peak memory bandwidth anyway so we don't lose anything
        # might be tempting to do chunks of 12 bytes, but that runs the risk of splitting across channels 
        # and also makes the data end up streaming slower (taking 2 cycles instead of 1)
        for x in range(0,(numBytes//8) + 1):
            ms = min(8,numBytes)
            self.memory.requestData(self,C,ms)
            self.memoryAccessBytes += ms
            numBytes -= ms
            # increment offset. If we move to a new channel, increment the channel number mod 8
            O += 8
            O %= 64
            if O < 8:
                C += 1
                C %= 8
        self.rowsToFetch = numRows # when this hits 0, the B prefetcher is done loading all data
    
    def cycle(self):
        if self.endFlag or not self.inputFlag:
            self.wastedCycles += 1
            return
        
        if self.currentRowSendingFinished and self.nextRowLoadingFinished:
            self.takeTimeStep()    
        elif self.currentRowSendingFinished and not self.nextRowLoadingFinished:
            self.memoryWastedCycles += 1
        

        # This is the part where we load a set of values into the multiplier array. Values are loaded 16 at a time
        # if our row pointer has exceeded the current length of the data, we don't do anything
        if not self.currentRowSendingFinished:
            end = self.currentRowPointer + 16 if self.currentRowPointer + 16 < len(self.currentRowData) else len(self.currentRowData)
            self.multiplierArray.loadBrow(self.currentRowData[self.currentRowPointer:end], self.currentRowIndices[self.currentRowPointer:end])
            self.currentRowPointer += 16
            # if we surpass the length of the current row after sending this data, we pop the row's row distance value and set a boolean
            # representing that the loading is finished and we can move on to the next row if the next set is in the buffer.
            if self.currentRowPointer >= len(self.currentRowData):
                self.rowDistanceHash[self.currentRowNumber].popleft()
                self.currentRowSendingFinished = True
    
        # This is the part where we load a row into buffer
        # MemoryFlag is false when prefetcher is waiting for memory, and true when it isnt
        if not self.nextRowLoadingFinished and self.rowsToFetch <= 0:                
            if self.nextRowPointer >= len(self.nextRowData):
                self.nextRowLoadingFinished = True
                return
            # Checks if the furthest row still has values to replace
            if self.prefetchBufferHash[self.furthestRowInPrefetcher] != 0:                                                
                # the cap on the number of cycles it'll take is either the amount of open space in the furthest row of the prefetcher
                # or, the amount of data actually needed to fill the row into the prefetcher.  We do a memory load of this size, 
                # and if we finish loading after this memory access, we set nextRowLoadingFinished to True. 
                # If we don't, then we come back and find the next lowest row in the buffer to replace.
                
                # ALSO: we will keep coming back and checking the next lowest row to replace because the furthest row in prefetcher could be the
                # current row, so we have to wait until its finished loading and gets popped before we indicate that it is after our current next row
                # and can start loading into its position.
                
                
                # Assumption made here: one element = one value/coordinate pair from the B matrix. 48 elements in each line.
                amountDataNeeded = min(self.prefetchBufferHash[self.furthestRowInPrefetcher] * 48, len(self.nextRowData) - self.nextRowPointer)
                
                self.nextRowPointer += amountDataNeeded # increment the row pointer by the amount of elements loaded into memory
                # Decrement the number of lines that are "eaten up" that were spilt
                amountDataNeeded = math.ceil(amountDataNeeded/48) # this is now in terms of the number of rows needed to be fetched
                self.prefetchBufferHash[self.furthestRowInPrefetcher] -= amountDataNeeded
                # Increment the number of lines that "replace" the previous row that was there
                if self.nextRowNumber in self.prefetchBufferHash:
                    self.prefetchBufferHash[self.nextRowNumber] += amountDataNeeded
                else:
                    self.prefetchBufferHash[self.nextRowNumber] = amountDataNeeded
                
                self.memoryAccess(amountDataNeeded, self.nextRowNumber, self.nextRowPointer) # pull indices and data into the buffer, this is counted in terms of prefetch buffer lines
            else:
                # this code runs when the fetchers need a new place to fill in the buffer, it selects a new furthest row, and deletes the current
                # furthest row in the prefetchBufferHash
                self.prefetchBufferHash.pop(self.furthestRowInPrefetcher)
                
                maximum = -1
                for key in self.prefetchBufferHash:
                    v = self.rowDistanceHash[key]
                    if not v:
                        maximum = key
                        break
                    if self.rowDistanceHash[maximum] == []:
                        break
                    if v[0] > self.rowDistanceHash[maximum][0]:
                        maximum = key
                
                self.furthestRowInPrefetcher = maximum
        elif self.rowsToFetch > 0:
            # reduce the number of rows to fetch by one for every 48 bytes that get fetched
            # since we don't use the bytes until the entire row is fetched, the order doesn't really matter anyway
            while self.freeBytes >= 48 * 12:
                self.rowsToFetch -= 1
                self.freeBytes -= 48 * 12
            return
                               

class MultiplierArray:
    def __init__(self) -> None:
        self.endFlag = False
        self.Avalue = 0
        self.AcompressedColumn = 0
        self.ARow = 0
        # we don't need the original column of A, but we do need a set of columns for B
        self.Bcols = []
        self.Brow = []
        
        # now, we need to make sure that multiplication will happen at the right cycle. Thus, all values are stored in "temp" values at first,
        # and then used next cycle
        self.Temp_Avalue = 0
        self.Temp_AcompressedColumn = 0
        self.Temp_ARow = 0
        # we don't need the original column of A, but we do need a set of columns for B
        self.Temp_Bcols = []
        self.Temp_Brow = []
        
        self.wastedCycles = 0

    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    # This runs every cycle that multiplication happens    
    def loadBrow(self,Brow,Bcols):
        self.Temp_Brow = Brow
        self.Temp_Bcols = Bcols
        #print(Brow, Bcols)
    
    def loadAvalue(self,AData):
        # no need for an EOF condition or anything here, if there isn't any input from B, there also isn't any input from A for a new row. When it starts back up again,
        # A will be sent along with the first row of B
        self.Temp_Avalue = AData[0]
        self.Temp_AcompressedColumn = AData[3]
        self.Temp_ARow = AData[1]
        #print(AData)
    
    # ensures that input is taken in next cycle, not the current one. Runs after all cycles are completed, is called by the outer function at the end of every loop.
    def loadInputs(self):
        self.Brow = self.Temp_Brow
        self.Bcols = self.Temp_Bcols
        self.Avalue = self.Temp_Avalue
        self.ARow = self.Temp_ARow
        self.AcompressedColumn = self.Temp_AcompressedColumn
        
        self.Temp_Bcols = []
        self.Temp_Brow = []
        
    
    def setMergeTree(self,MT):
        self.MergeTree = MT
    
    def setScheduler(self,S:"Scheduler"):
        self.scheduler = S
    
    
    def cycle(self):
        # this acts as the "input flag" that we have on other parts
        if self.scheduler.merging or self.endFlag:
            self.wastedCycles += 1
            return

        if self.Brow:
            output = []
            for x in range(len(self.Brow)):
                output.append((self.Brow[x]*self.Avalue,self.ARow,self.Bcols[x])) # send value, row, column in COO format
            self.Brow = [] #set it to empty, we need new input in order to multiply again.
            self.Bcols = []
            self.MergeTree.input(output, self.AcompressedColumn) # AcompressedColumn indicates which input buffer this needs to go into
        else:
            self.scheduler.checkEOF()
            
             
class MergeTree:
    # 1) This class contains the huffman tree scheduler using a minheap (to choose which matrices to merge)
    # When the merge tree is empty, since the prio queue can be added to easily (summing up all components), we run the huffman tree scheduler and notify
    # the partial matrix fetcher as well as the matrix A loader
    # 2) This class contains the merge tree (starting from 64, 32, 16,...)
    
    # Each fifo has to have 48 values for the following reasons:
    # 1) each merger in the merge tree can do a single 16-value merge each cycle. This means that the cap of the fifo has to be the highest
    # number achievable so that a merger does not get underutilized and be unable to keep up with the throughput of the overall tree
    # 2) a merger cannot react to values pulled in the same cycle. It needs to wait for the next cycle to fill up the values taken out
    # by its upper level merger in the previous cycle.
    def __init__(self, mergeSize) -> None:
        self.merging = False
        self.mergeSize = mergeSize
        self.loaded = [True] + [False] * 6
        self.inputFifos = [[deque([]),False] for x in range(64)]
        self.L1 = [[deque([]),True] for _ in range(32)]
        self.L2 = [[deque([]),True] for _ in range(16)]
        self.L3 = [[deque([]),True] for _ in range(8)]
        self.L4 = [[deque([]),True] for _ in range(4)]
        self.L5 = [[deque([]),True] for _ in range(2)]
        self.outputFifo = [[deque([]),True]]
        self.fifos = [self.inputFifos, self.L1, self.L2, self.L3, self.L4, self.L5, self.outputFifo]
        self.out = []
        self.scheduler = None
        self.totalCount = 0
        self.endFlag = False
        
        # The partial matrix storage is also implemented in this class. Yes, it isn't like that in the actual paper, but
        # it doesn't really need its own class since the only thing it interacts with is the Merge Tree anyway
        
        self.loadingPartials = False
        self.partialMatrices = []
        self.memoryLatency = 0
        self.wastedCycles = 0
        self.idleCycles = 0
        self.totalMergeCycles = 0
        self.totalEmptyFifos = 0
    
    def reset(self):
        self.merging = False
        self.loaded = [True] + [False] * 6
        self.inputFifos = [[deque([]),False] for x in range(64)]
        self.L1 = [[deque([]),True] for _ in range(32)]
        self.L2 = [[deque([]),True] for _ in range(16)]
        self.L3 = [[deque([]),True] for _ in range(8)]
        self.L4 = [[deque([]),True] for _ in range(4)]
        self.L5 = [[deque([]),True] for _ in range(2)]
        self.outputFifo = [[deque([]),True]]
        self.fifos = [self.inputFifos, self.L1, self.L2, self.L3, self.L4, self.L5, self.outputFifo]
        self.out = []
    
    def input(self, Mvals, cc):
        self.inputFifos[cc][0].extend(Mvals)
        
    def setScheduler(self, S: "Scheduler"):
        self.scheduler = S
        
    def startMerging(self):
        self.merging = True
        
    def inputPartials(self, num, matrices, newMatrixSize):
        self.loadingPartials = True
        self.memoryLatency = 0
        for x in matrices:
            self.memoryLatency += 0
        self.totalCount = newMatrixSize
        print(len(matrices))        
        for i,x in enumerate(matrices):
            self.inputFifos[i + num][0].extend(self.partialMatrices[x])

    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def cycle(self):
        if self.endFlag:
            self.wastedCycles += 1
            return
        if self.loadingPartials:
            if self.memoryLatency <= 0:
                self.loadingPartials = False
            else:
                self.memoryLatency -= 128
        elif not self.merging:
            # not loading partials or merging, so its a wasted cycle
            self.wastedCycles += 1
            return
        if not self.merging:
            # this can't be merged with the EndFlag above, because we need to load partials while not merging
            # loading partials and not merging, idle cycle
            self.idleCycles += 1
            return

        
        # first, we choose which queue we will fill up in each row
        emptyFifos = 0
        chosenQueues = []

        for i in range(1,len(self.fifos)):
            # Checks whether a non fully-loaded layer can be registered as fully loaded again
            if not self.loaded[i]:
                newVal = True
                for queue in self.fifos[i]:
                    if len(queue[0]) < self.mergeSize and queue[1]:
                        newVal = False
                    elif len(queue[0]) == 0 and not queue[1]:
                        emptyFifos += 1
                self.loaded[i] = newVal

            if self.loaded[i-1]:
                minimum = 10000000 # unrealistic, impossible number
                queueToFill = -1
                for num,queue in enumerate(self.fifos[i]):                    
                    # conditions: length of queue is shorter than current lowest, the queue has not "dried up", 
                    # and its two lower values are either dried up or have more than 16 values inside them.
                    if len(queue[0]) < minimum and queue[1] and (len(self.fifos[i-1][num*2][0]) > self.mergeSize or not self.fifos[i-1][num*2][1]) and (len(self.fifos[i-1][num*2 + 1][0]) > self.mergeSize or not self.fifos[i-1][num*2 + 1][1]):
                        minimum = len(queue)
                        queueToFill = num
                    elif len(queue[0]) == 0 and not queue[1]:
                        # this means that in this cycle the queue is empty, and will not be used again.
                        # we don't consider queues that could be used again for this, because they could only be empty for this cycle
                        # and are very much an important part of the computation (they are being utilized in some way)
                        emptyFifos += 1
                chosenQueues.append(queueToFill) # the queue we will be filling using the ith comparator array is queueToFill
            else:
                # if the subqueues have not been loaded yet:
                chosenQueues.append(None) 
        self.totalMergeCycles += 1
        self.totalEmptyFifos += emptyFifos/64
        # in this iteration, a few things can happen:
        # 1) the prior layer is not full enough to begin merging (chosenQueues has None)
        # 2) the merger was unable to choose a queue to fill (chosenQueues has -1). In this case, we don't do anything
        # 3) the merger chose a queue to fill. In this case, we pull values from the previous layer to fill the queue, and then check 
        #    the previous layer's fifos to see if we need to check this queue off as having streamed everything it could 
        #    (the queues that feed into it are empty and also have streamed everything)
        for i,q in enumerate(chosenQueues):
            if q == None:
                break #none of the other queues will have enough either
            elif q == -1:
                # we don't merge anything for this layer. either the previous layer doesn't have enough values to make safe merges happen
                # or, the entire layer is finished computing (in which case nothing needs to happen)
                continue 
            else:
                # this part is not done using the fancy merger, I abstract it out and just use a pointer-based merger instead
                values1 = self.fifos[i][q*2][0] #remember, i goes from 0, so we don't need to subtract 1 to get the previous layer
                values2 = self.fifos[i][q*2 + 1][0]
                count = 0
            
                while values1 and values2 and count < self.mergeSize:
                    c1 = values1[0]
                    c2 = values2[0]
                    if c1[1] < c2[1]: # if the row value of 1 is greater
                        self.fifos[i+1][q][0].append(values1.popleft())
                    elif c1[1] > c2[1]: # if the row value of 2 is greater
                        self.fifos[i+1][q][0].append(values2.popleft())
                    elif c1[2] < c2[2]: # if the coulumn value of 1 is greater
                        self.fifos[i+1][q][0].append(values1.popleft())
                    elif c1[2] > c2[2]: # if the column value of 2 is greater
                        self.fifos[i+1][q][0].append(values2.popleft())
                    else: # they must be equal
                        self.fifos[i+1][q][0].append((c1[0]+c2[0],c1[1],c1[2]))
                        values1.popleft()
                        values2.popleft()
                    count += 1
                
                # in case one row ran out (meaning all previous fifos under it are also empty!)
                while values1 and count < self.mergeSize:
                    self.fifos[i+1][q][0].append(values1.popleft())
                    count += 1
                while values2 and count < self.mergeSize:
                    self.fifos[i+1][q][0].append(values2.popleft())
                    count += 1
                
                # this is where we check whether we shut down this FIFO or not
                # The conditions are if its two children have verified that all of their children are empty, 
                # and if its two children are themselves empty. We verify this by the boolean attached to each fifo.
                # note that when initially filling in a new layer, if we aren't merging at full capacity, 
                # all the empty inner nodes will be shut down as they try to load from the previous layer.
                if not self.fifos[i][q*2][0] and not self.fifos[i][q*2 + 1][0] and not self.fifos[i][q*2][1] and not self.fifos[i][q*2 + 1][1]:
                    self.fifos[i+1][q][1] = False

        if len(self.fifos[-1][0][0]) >= self.mergeSize or not self.fifos[-1][0][1]:
            for x in range(min(self.mergeSize,len(self.fifos[-1][0][0]))):
                self.out.append(self.fifos[-1][0][0].popleft()) # pop 16 values out of the root fifo if possible.
        if not self.fifos[-1][0][1] and not self.fifos[-1][0][0]: # End of merging
            self.merging = False
            self.partialMatrices.append(self.out) # add output matrix to the partial matrix storage
            self.reset()
            self.scheduler.addMatrix(len(self.partialMatrices)-1,self.totalCount) # store new matrix in scheduler
            self.scheduler.schedule()
