import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from collections import deque
import math
import sys



class Scheduler:
    # The scheduler is where the Huffman Tree is located. When the MergeTree finishes its merging, it calculates the next set of partial matrices to merge
    # it then goes to the PartialMatrixFetcher to pull any partial matrices that are small enough into the mergeTree, and to the MatrixAFetcher to pull 
    # a set of compressed columns (each of which becomes a partial matrix to be sent to the MergeTree). 
    def __init__(self, MAF: "MatrixAFetcher", MT:"MergeTree", DLB: "DistanceListBuilder", MBP: "MatrixBPrefetcher") -> None:
        self.MatrixAFetcher = MAF
        self.MergeTree = MT
        self.DistanceListBuilder = DLB
        self.MatrixBPrefetcher = MBP
        self.MatrixARows = MAF.rowLengths.copy()
        self.StoredMatrices = []
        self.maxSchedule = ((len(self.MatrixARows)-2) % (64-1)) + 2
        self.merging = False

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
        
        print(self.StoredMatrices)
        self.maxSchedule = 64
        self.MatrixAFetcher.inputInstructions(numACols)
        self.MergeTree.inputPartials(numACols,matrices,PMlength)
        self.merging = False
        
        
    # Called by the MultiplierArray when it has no input
    def checkEOF(self):
        # print(self.merging, self.MatrixAFetcher.inputFlag, self.DistanceListBuilder.lookaheadFIFO, self.MatrixBPrefetcher.inputFlag)
        # if not merging, and no input remaining:
        if not self.merging and not self.MatrixAFetcher.inputFlag and not self.DistanceListBuilder.lookaheadFIFO and not self.MatrixBPrefetcher.inputFlag:
            self.merging = True
            # indicates to the merge Tree that it should begin merging
            
            self.MergeTree.startMerging()
        
        
        
            
        
        

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
        self.rowLengths = []
        
        # preprocessing to make sure that we have a list of the number of values in each compressed column - important for the huffman tree scheduler
        for x in range(len(A.indptr)-1):
            l = A.indptr[x+1] - A.indptr[x]
            if len(self.rowLengths) < l:
                self.rowLengths += ([0] * (l - len(self.rowLengths)))
            for y in range(l):
                self.rowLengths[y] += 1
        # The currentCondensedColumn should be equal to the furthest condensed column at first
        self.currentCondensedColumn = len(self.rowLengths)
        self.currentNumColumns = 0
        self.currentRow = 0

    def setDLB(self,DLB):
        self.distanceListBuilder = DLB
    
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

        # Make use of current condensed column number, and set currentCondensedColumn (the minimum number of columns needed in a row to have values we want for this round)
        # currentNumColumns allows us to specify the range of columns we want to pull from memory (condensed columns [currendCondensedColumns:currentCondensedColumns + currentNumColumns])
        # This number is specified by the huffman tree builder, making use of the number of values in each column within a minheap
    
    def cycle(self):
        #Every cycle, if we arent at EOF:
        # 1) Look through a row, check if it has have columns within the condensed column.
        # 2) If it does, load them into the FIFO and increment
        # 3) If it doesn't, increment.
        if self.endFlag or not self.inputFlag:
            return
        
        if self.currentRow >= len(self.A.indptr) -1:
            # We have gone through all the rows, and all valid A values are sent. Send none so that the Multiplier can signal that loading is finished to the MergeTree
            self.inputFlag = False
            return
        
        if self.A.indptr[self.currentRow+1] - self.A.indptr[self.currentRow] > self.currentCondensedColumn:
            i = self.A.indptr[self.currentRow] + self.currentCondensedColumn
            #print(self.A.indptr[self.currentRow+1] - self.A.indptr[self.currentRow] - self.currentCondensedColumn)
            if self.A.indptr[self.currentRow+1] - self.A.indptr[self.currentRow] - self.currentCondensedColumn < self.currentNumColumns: 
                end = self.A.indptr[self.currentRow+1]
            else:
                end = self.currentNumColumns + self.A.indptr[self.currentRow] + self.currentCondensedColumn
            #print(i, end)
            for x in range(i,end):
                #print(i,x,x-i)
                self.distanceListBuilder.input((self.A.data[x],self.currentRow,self.A.indices[x],x-i)) # data, row, column, condensed column number are sent to DLB buffer. 
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

    def cycle(self):
        # This part of the accelerator really doesn't need a cycle-by-cycle thing, it just sends data to the matrix B prefetcher and 
        # multiplier array when asked. No point in modeling something that doesn't change anything
        pass
        
    
class MatrixBPrefetcher:
    # Data is sent to multiplier in chunks of 16 per cycle, once a row is sent, the new A-value is sent to the multiplier, kicking out the old one.
    # Then, the matrixBPrefetcher sends the values for the new row to the multiplier
    
    
    # An important note about this part of the accelerator is that the 16 data fetchers work in parallel, and they work to load new buffer lines as the preloader sends data to the 
    # multiplier array every cycle
    
    # Note: 
    # Option 1) when it says multiple fetchers, does it mean that individual rows are fetched in parallel, hiding the latency for each row?
    # Option 2) Maybe, instead of accessing the DRAM in parallel to load multiple rows at once, it loads in each row individually, and moves onto the next one before the first one is done multiplying
    # This way, there is only one "fetcher", but it loads in parallel?
    # Option 3) Each time the prefetcher is done sending a row to the multiplier, the buffer is "reformatted", and each fetcher is assigned a unique row to start/continue loading in
    # This seems needlessly complicated, and is not necessarily indicative of what is shown in Figure 9. 
    
    # As it stands, in Figure 9, no "prefetching" is going on. At every time step when a value is needed, it is brought in from memory
    
    
    
    
    """
    Alright nebil - here's the plan:
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
        self.memoryFlag = True
        self.distanceFIFO = deque([])
        self.prefetchBufferHash = {-1: 1024} # initially, every single set in the row is empty (assigned to -1), there are 1024 lines
        self.rowDistanceHash = {-1: deque([0])} # a queue of row distances in order assigned to this row
        self.multiplierArray = None
        
        
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
        self.memLatency = 0
        self.memoryAccessBytes = 0


    def takeTimeStep(self):
        # set this true, if we see that there is input in the future, we set it back to false
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
        
        
    def memoryAccess(self, numRows):
        if numRows == 0:
            return
        # we make use of 16 fetchers, so assuming each fetcher takes about the same time to load a row, ceil(numbytes/16)
        # Each row is 48 elements, and each element is 12 bytes
        self.memoryAccessBytes += numRows * 48 * 12 
        
        # The number of rows that need to be fetched by each fetcher:
        numRowsFetched = math.ceil(numRows/16)
        self.memLatency = 70 + 5 + math.ceil((numRowsFetched * 48 * 12)/128) # Memory Row latency of 70 cycles, Column latency of 5 cycles + 128 bytes/cycle, given 128gb/s and 1GhZ clock speed
        self.memoryFlag = True
    
    def cycle(self):
        if self.endFlag or not self.inputFlag:
            return
        
        if self.currentRowSendingFinished and self.nextRowLoadingFinished:
            self.takeTimeStep()

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
        if not self.nextRowLoadingFinished and not self.memoryFlag:                
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
                
                self.memoryAccess(amountDataNeeded) # pull indices and data into the buffer, this is counted in terms of prefetch buffer lines
            else:
                # this code runs when the fetchers need a new place to fill in the buffer, it selects a new furthest row, and deletes the current
                # furthest row in the prefetchBufferHash
                self.prefetchBufferHash.pop(self.furthestRowInPrefetcher)
                
                maximum = -1
                for key in self.prefetchBufferHash:
                    v = self.rowDistanceHash[key]
                    if v == []:
                        maximum = key
                        break
                    if self.rowDistanceHash[maximum] == []:
                        break
                    if v[0] > self.rowDistanceHash[maximum][0]:
                        maximum = key
                
                self.furthestRowInPrefetcher = maximum
        elif self.memoryFlag:
            self.memoryWastedCycles += 1
            self.memLatency -= 1
            if self.memLatency <= 0:
                self.memoryFlag = False
            return
                               

class MultiplierArray:
    def __init__(self) -> None:
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
        if self.scheduler.merging:
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
    def __init__(self) -> None:
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
        self.scheduler = None
        self.totalCount = 0
        
        # The partial matrix storage is also implemented in this class. Yes, it isn't like that in the actual paper, but
        # it doesn't really need its own class since the only thing it interacts with is the Merge Tree anyway
        
        self.loadingPartials = False
        self.partialMatrices = []
        self.memoryLatency = 0
    
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
        self.memoryLatency = 0 #this is where the memory latency calculations would be done,
        self.totalCount = newMatrixSize
        print(matrices)        
        for i,x in enumerate(matrices):
            self.inputFifos[i + num][0].extend(self.partialMatrices[x])
    
    def cycle(self):
        if self.loadingPartials:
            if self.memoryLatency <= 0:
                self.loadingPartials = False
            else:
                return
        if not self.merging:
            return
        
        # first, we choose which queue we will fill up in each row
        chosenQueues = []

        for i in range(1,len(self.fifos)):
            # Checks whether a non fully-loaded layer can be registered as fully loaded again
            if not self.loaded[i]:
                newVal = True
                for queue in self.fifos[i]:
                    if len(queue[0]) < 16 and queue[1]:
                        newVal = False
                        break
                self.loaded[i] = newVal

            if self.loaded[i-1]:
                minimum = 10000000 # unrealistic, impossible number
                queueToFill = -1
                for num,queue in enumerate(self.fifos[i]):                    
                    
                    # conditions: length of queue is shorter than current lowest, the queue has not "dried up", 
                    # and its two lower values are either dried up or have more than 16 values inside them.
                    if len(queue[0]) < minimum and queue[1] and (len(self.fifos[i-1][num*2][0]) > 16 or not self.fifos[i-1][num*2][1]) and (len(self.fifos[i-1][num*2 + 1][0]) > 16 or not self.fifos[i-1][num*2 + 1][1]):
                        minimum = len(queue)
                        queueToFill = num
                chosenQueues.append(queueToFill) # the queue we will be filling using the ith comparator array is queueToFill
            else:
                # if the subqueues have not been loaded yet:
                chosenQueues.append(None) 
        
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
            
                while values1 and values2 and count < 16:
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
                while values1 and count < 16:
                    self.fifos[i+1][q][0].append(values1.popleft())
                    count += 1
                while values2 and count < 16:
                    self.fifos[i+1][q][0].append(values2.popleft())
                    count += 1
                
                # this is where we check whether we shut down this FIFO or not
                # The conditions are if its two children have verified that all of their children are empty, 
                # and if its two children are themselves empty. We verify this by the boolean attached to each fifo.
                # note that when initially filling in a new layer, if we aren't merging at full capacity, 
                # all the empty inner nodes will be shut down as they try to load from the previous layer.
                if not self.fifos[i][q*2][0] and not self.fifos[i][q*2 + 1][0] and not self.fifos[i][q*2][1] and not self.fifos[i][q*2 + 1][1]:
                    self.fifos[i+1][q][1] = False

        if len(self.fifos[-1][0][0]) >= 16 or not self.fifos[-1][0][1]:
            for x in range(min(16,len(self.fifos[-1][0][0]))):
                self.out.append(self.fifos[-1][0][0].popleft()) # pop 16 values out of the root fifo if possible.
        if not self.fifos[-1][0][1] and not self.fifos[-1][0][0]: # End of merging
            self.merging = False
            self.partialMatrices.append(self.out) # add output matrix to the partial matrix storage
            self.reset()
            self.scheduler.addMatrix(len(self.partialMatrices)-1,self.totalCount) # store new matrix in scheduler
            self.scheduler.schedule()
        

"""
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from collections import deque
import math

global total_memory_usage
global truepartial

truepartial = []

total_memory_usage = 0

#Matrix Condensing
class Afetcher:
    def __init__(self, A) -> None:
        global total_memory_usage 

        self.condenseSize = 63
        self.A = csr_matrix(A)
        data = self.A.data
        indices = self.A.indices
        self.originalcols = np.copy(self.A.indices)
        indptr = self.A.indptr
        
        
        # the entire A matrix is fetched from memory column by column. I can just calculate how much data is pulled all at once.
        total_memory_usage += len(data) + len(indptr) + len(indices)

        
        # rearrange the order of the indices so that it is in increasing sorted order. Might be unnecessary but it will break otherwise
        for i in range(0,len(indptr)-1):
            data[indptr[i]:indptr[i+1]] = data[indptr[i] + np.argsort(indices[indptr[i]:indptr[i+1]])]
            self.originalcols[indptr[i]:indptr[i+1]] = self.originalcols[indptr[i] + np.argsort(indices[indptr[i]:indptr[i+1]])]
            indices[indptr[i]:indptr[i+1]] = np.arange(0,indptr[i+1]-indptr[i])
        
        self.A = csr_matrix((data,indices,indptr))
        self.colPointer = self.A.shape[1]        
        self.rowPointer = 1
    
    def __str__(self) -> str:
        return str(self.A.toarray()) + "\n column pointer at: " + str(self.colPointer)
    

    def resetRows(self):
        self.rowPointer = 1
        self.colPointer -= self.condenseSize
    
    #returns a set of values that are in the same row of the current compressed column.
    def nextRowSet(self):
        
        if self.rowPointer > self.A.shape[0]:
            return None
            
        minval = 0 if self.colPointer - self.condenseSize <= 0 else self.colPointer-self.condenseSize
        rp0 = self.A.indptr[self.rowPointer-1]
        rp1 = self.A.indptr[self.rowPointer]

        mask = np.in1d(self.A.indices[rp0: rp1],np.arange(minval, self.colPointer))
        
        while not mask.any() and self.colPointer > 0:
            if self.rowPointer >= self.A.shape[0]:
                return None
            else:
                self.rowPointer += 1
            if self.colPointer < 0:
                return None
            
            minval = 0 if self.colPointer - self.condenseSize <= 0 else self.colPointer-self.condenseSize

                
            rp0 = self.A.indptr[self.rowPointer-1]
            rp1 = self.A.indptr[self.rowPointer]
            mask = np.in1d(self.A.indices[rp0: rp1],np.arange(minval, self.colPointer))
        
        compressedcols = list(map(lambda x: x%self.condenseSize, self.A.indices[rp0:rp1][mask]))
        retval = [self.rowPointer-1, compressedcols,  self.originalcols[rp0:rp1][mask],self.A.data[rp0:rp1][mask]]
        
        self.rowPointer += 1
        return retval

# Row Prefetcher
class Bprefetcher:
    def __init__(self,B) -> None:
        self.B = csr_matrix(B)
    def fetch(self, row):
        # Note: I tried to implement the B prefetcher, but it got way too complex so i decided to just abstract it out.
        # Since the row prefetcher was stated to achieve a 62% hit rate, I'm just going to count the amount of data pulled from memory and then 
        # multiply it by 0.38 to save some headache.
        global total_memory_usage
  
        begin = self.B.indptr[row]
        end = self.B.indptr[row+1]
        data = [self.B.data[begin:end],self.B.indices[begin:end]]
        total_memory_usage += math.ceil(len(data) * 0.38)
        return data


class Merger:
    def __init__(self, fifolist,nextMerger) -> None:
        self.fifolist = fifolist
        self.nextMerger = nextMerger
    def __bool__(self):
        for x in self.fifolist:
            if bool(x):
                return True
        return False

    def __str__(self) -> str:
        ret = ""
        for x in self.fifolist:
            if bool(x):
                ret += str(x)        
        return ret
        
    def MergeAndPush(self, ind):
        if self.nextMerger == None:
            global truepartial
            truepartial += list(self.fifolist[0])
            self.fifolist[0] = deque()
            return
    
        i1 = self.fifolist[ind * 2]
        i2 = self.fifolist[ind * 2 + 1]
        c1 = 0
        c2 = 0
        while bool(i1) and bool(i2) and c1 < 16 and c2 < 16:
            if i1[0][0] < i2[0][0]:
                self.nextMerger.fifolist[ind].append(i1.popleft())
                c1 += 1
            elif i1[0][0] > i2[0][0]:
                self.nextMerger.fifolist[ind].append(i2.popleft())
                c2 += 1
            elif i1[0][1] < i2[0][1]:
                self.nextMerger.fifolist[ind].append(i1.popleft())
                c1 += 1
            elif i1[0][1] > i2[0][1]:
                self.nextMerger.fifolist[ind].append(i2.popleft())
                c2 += 1
            else:
        
                self.nextMerger.fifolist[ind].append([i1[0][0],i1[0][1],i1.popleft()[2] + i2.popleft()[2]])
                c1 += 1
                c2 += 1
        if c1 == 16 or not bool(i1):
            while bool(i2) and c2 < 16:
                self.nextMerger.fifolist[ind].append(i2.popleft())
                c2 += 1
        elif c2 == 16 or not bool(i2):
            while bool(i1) and c1 < 16:
                self.nextMerger.fifolist[ind].append(i1.popleft())
                c1 += 1
        
            


class MultiplyAndMerge:
    def __init__(self, A, B) -> None:
        global truepartial
        self.aFetcher = Afetcher(A)
        self.bPrefetcher = Bprefetcher(B)
        self.partials = [deque() for _ in range(64)]
        for x in truepartial:
            self.partials[63].append(x)
        self.endflag = False
        self.numcycles = 0
        
    def MultiplyColumn(self):
        rs = self.aFetcher.nextRowSet()
        while rs != None:
            for x in range(len(rs[3])):
                brow = np.array(self.bPrefetcher.fetch(rs[2][x]))
                datapoints = brow[0] * rs[3][x]
                self.numcycles += len(brow) # vectorized multiply across all of B
                for val in range(len(datapoints)):
                    self.partials[rs[1][x]].append([rs[0], brow[1][val],datapoints[val]])
            rs = self.aFetcher.nextRowSet()
        self.aFetcher.resetRows()
        
    
    def MergePartials(self):
        self.endflag = False
        out = Merger([deque()],None)
        l6m = Merger([deque() for _ in range(2)],out)
        l5m = Merger([deque() for _ in range(4)],l6m)
        l4m = Merger([deque() for _ in range(8)],l5m)
        l3m = Merger([deque() for _ in range(16)],l4m)
        l2m = Merger([deque() for _ in range(32)],l3m)
        l1m = Merger(self.partials,l2m)
        
        mergerArray = [out,l6m,l5m,l4m,l3m,l2m,l1m]
        
        
        while bool(out.fifolist[0]) or not self.endflag:
            self.numcycles += 1
            out.MergeAndPush(0)
            for p in range(6):
                for val in range(2**(p)):
                    c = mergerArray[p+1]
                    n = mergerArray[p]
                    # if the fifo in the next level is empty:
                    if not bool(n.fifolist[val]) and (bool(c.fifolist[val*2]) or bool(c.fifolist[val*2+1])):
                        c.MergeAndPush(val)
                        break
            # if the out fifo is empty, we check each layer to see if the entire thing is empty
            if not bool(out.fifolist[0]):
                for x in mergerArray:
                    if bool(x):
                        break
                else:
                    self.endflag = True




def run_sparch(matrix1,matrix2):
    
    global total_memory_usage

    merger = MultiplyAndMerge(matrix1,matrix2)
    merger.MultiplyColumn()
    merger.MergePartials()
    
    
    data = []
    i = []
    j = []
    
    global truepartial
    for x in truepartial:
        i.append(x[0])
        j.append(x[1])
        data.append(x[2])
    
    truepartial = []

    output = coo_matrix((data, (i,j)), shape= (matrix1.shape[0], matrix2.shape[1]))

    print("number of cycles:", merger.numcycles )
    print("Data Pulled from DRAM:", total_memory_usage)
    if matrix1.size >= 10000 or matrix2.size >= 10000:
        print("matrices too large to verify")
    #else:
    #    trueval = matrix1 @ matrix2
    #    print("Verify that sparse multiplication is correct: ", np.allclose(output.toarray(),trueval,rtol=0.000001))
    
    total_memory_usage = 0
    
"""