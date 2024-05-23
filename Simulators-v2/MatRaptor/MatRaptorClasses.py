import numpy as np
from collections import deque
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

import scipy.io as sio
import scipy.spatial

# Step 1: Convert CSR to C2SR - C2SR basically just makes sure that all the data a PE reads at a given time is useful.
# Notes on how I implemented C2SR:
# 1) the row lengths and row pointers point to CHANNELS where data is kept, and the index location in that channel. This is honestly to avoid
# having to split up a single row among many different parts of the data and column id arrays. Instead of using one list, I just use the number
# of lists equivalent to the number of channels and then make a list for each one since they are vectorized anyway it doesn't change anything
# and makes the code less confusing later on.
def csr_to_c2sr(data, indices, indptr, num_channels):
    # Create arrays to store C2SR format data
    c2sr_values = [[] for _ in range(num_channels)]
    c2sr_column_ids = [[] for _ in range(num_channels)]
    c2sr_row_length = []
    c2sr_row_pointer = []

    # Loop through each row in CSR format
    for i in range(len(indptr) - 1):
        # Get the start and end indices for the current row in CSR format
        start_idx = indptr[i]
        end_idx = indptr[i + 1]

        # Get the channel assigned to the current row
        channel = i % num_channels
        
        # Append the row length and row pointer for the current row in C2SR format
        c2sr_row_length.append((channel, end_idx - start_idx))
        c2sr_row_pointer.append((channel, len(c2sr_values[channel])))

        # Loop through the non-zero elements of the current row
        for j in range(start_idx, end_idx):
            # Append the value and column id for the current element in C2SR format
            c2sr_values[channel].append(data[j])
            c2sr_column_ids[channel].append(indices[j])
        
        
        

    return c2sr_values, c2sr_column_ids, c2sr_row_length, c2sr_row_pointer

class PE:
    # Notes about PEs:
    # 1 - Input comes in the form of (a,b, i, j) pairs, where "k" is the missing coordinate (col of a, row of b)
    # Each time the "k" changes, the PE needs to find the smallest possible queue in order to merge the next set together. 
    # Because k is not passed into the PE, this switch will come in the form of a flag being set on the PE by the SpBL
    # 2 - The PE and SpBL will have buffers for their input.
    # 3 - Each PE will have a number of Queues, and two sets of these queues. The PE will switch between filling up one set of queues and merging/streaming the other to DRAM (as a complete row) when i changes
    # 4 - Look into queue overflow.
    # 5 - Basically all computation here happens in the PEs, need to be careful with the multiple sets of PEs, and properly using the adder tree.
    
    def __init__(self, numQueues, PENumber, mode) -> None:
        self.inputBuffer = deque([])
        self.outputBuffer = deque([])
        self.mode = mode
        
        
        if mode != "SmallMerge":
            self.inputQueues = []
            self.outputQueues = []
            for x in range(0,numQueues):
                self.inputQueues.append(deque([],maxlen=20000))
                self.outputQueues.append(deque([],maxlen=20000))
            self.inputQueuesLengths = [0] * numQueues # stores the length of each queue so the shortest one can be found efficiently.

            self.helperQueueNumber = numQueues-1
            self.currentQueueNumber = 0
            self.tempI = None # this is used to switch I values while keeping track of previous I value so it can be shared to part 2
            self.currentI = -1
            self.prevI = -1
            self.currentK = -1
            self.PENumber = PENumber
            self.numQueues = numQueues
            self.endFlag = False # indicates end of processing
            self.endFlag1 = False # indicates end of processing for part I
            self.endFlag2 = False # indicates the next swap is the end of processing for part II
            self.outputEmptyFlag = False # indicates that the outgoing queues are empty and ready to be swapped
            self.inputFinishedFlag = False # indicates that i (row of a) has been incremented, and the queue set is ready to be swapped
            self.inputFlag = False # indicates that the input queue is empty.
            self.numWastedCycles = 0
            self.partIWastedCycles = 0
            self.partIIWastedCycles = 0
        else:  
            self.inputQueues = [deque([]), deque([])]
            self.outputQueues = [deque([]), deque([])]
            self.helperQueueNumber = 0
            self.currentQueueNumber = 1
            self.outputQueueNumber = 0
            self.currentK = -1
            self.currentI = -1
            self.tempI = None # this is used to switch I values while keeping track of previous I value so it can be shared to part 2
            self.currentI = -1
            self.prevI = -1
            self.currentK = -1
            self.PENumber = PENumber
            self.numQueues = numQueues
            self.endFlag = False # indicates end of processing
            self.endFlag1 = False # indicates end of processing for part I
            self.endFlag2 = False # indicates the next swap is the end of processing for part II
            self.outputEmptyFlag = False # indicates that the outgoing queues are empty and ready to be swapped
            self.inputFinishedFlag = False # indicates that i (row of a) has been incremented, and the queue set is ready to be swapped
            self.inputFlag = False # indicates that the input queue is empty.
            self.numWastedCycles = 0
            self.partIWastedCycles = 0
            self.partIIWastedCycles = 0
        
        
        
    def __str__(self) -> str:
        s = ""
        s += "inputQueues: "
        for i,x in enumerate(self.inputQueues):
            s += str(i) + ": "
            s += str(x)
            s += ", "
            
        s += "\n outputQueues: "
        for i,x in enumerate(self.outputQueues):
            s += str(i) + ": "
            s += str(x)
            s += ", "
        
        s += "\n inputQueueLengths: "
        for i,x in enumerate(self.inputQueuesLengths):
            s += str(i) + ": "
            s += str(x)
            s += ", " 
        
        s += "\n inputFlag: " + str(self.inputFlag)
        s += "\n inputBuffer: " + str(self.inputBuffer)
        s += "\n prevI: " + str(self.prevI)
        s += "\n currentI: " + str(self.currentI)
        s += "\n currentQN: " + str(self.currentQueueNumber)
        s += "\n helperQN: " + str(self.helperQueueNumber)


        
        return s
    
    def input(self,a,b,i,j, k):
        self.inputBuffer.append((a,b,i,j, k))
        self.inputFlag = True
        
    def cycle(self):
        if not self.inputFlag or self.endFlag:
            # No Computation to Be Done.
            self.numWastedCycles += 1
            return
       

        # BEGIN SMALLMERGE
        if self.mode == "SmallMerge":
            #print(self.outputEmptyFlag, self.inputFinishedFlag)
            if self.outputEmptyFlag and self.inputFinishedFlag:

                #print("SWAP")
                self.outputEmptyFlag = False
                self.inputFinishedFlag = False
                
                if self.endFlag2:
                    self.endFlag = True
                    return
                if self.endFlag1: # If we are swapping and there is no more input, we need to swap one last time so the last numbers can be cleaned up
                    #print("END2")
                    self.endFlag2 = True
                    self.inputFinishedFlag = True
                    
                #swap queue sets, set appropriate values, reset flags
                temp = self.inputQueues
                self.inputQueues = self.outputQueues
                self.outputQueues = temp
                
                self.outputQueueNumber = self.currentQueueNumber

                self.helperQueueNumber = 0
                self.currentQueueNumber = 1
                
                self.prevI = self.currentI 
                self.currentI = self.tempI #remember, tempI is set when the inputFinishedFlag is set, with the new value of I
                
                # NOTE: at this point, self.inputQueues should be full of empty queues, and all of its helper variables should be reset
                # On the other hand, self.outputQueues should be full of values ready to be merged.
                
 
            # Part II
            if not self.outputEmptyFlag:
                
                # Check if prevI is valid!
                if self.prevI == -1:
                    self.outputEmptyFlag = True
                else:
                    if self.outputQueues[self.outputQueueNumber]:
                        v = self.outputQueues[self.outputQueueNumber].popleft()
                        self.outputBuffer.append((v[0],self.prevI,v[1]))
                    else:
                        self.outputEmptyFlag = True
            else:
                self.partIIWastedCycles += 1

            

            # Part I
            if not self.inputFinishedFlag or not self.endFlag1:
                if not self.inputBuffer:
                    self.inputFlag = False # If input buffer is empty, do not run
                    return
                
                a, b, i, j, k = self.inputBuffer.popleft()
                if a == None or k != self.currentK or i != self.currentI:
                    l = len(self.inputQueues[self.currentQueueNumber])
                    for x in range(l):
                        self.inputQueues[self.helperQueueNumber].append(self.inputQueues[self.currentQueueNumber].popleft())
                    temp = self.currentQueueNumber
                    self.currentQueueNumber = self.helperQueueNumber
                    self.helperQueueNumber = temp
                    self.currentK = k
                
                if a == None:
                    self.inputFlag = True
                    self.inputFinishedFlag = True
                    self.endFlag1 = True
                    return # End of Processing

                if i != self.currentI:
                    # If I has been changed, then set the inputFinishedFlag and change currentI, and skip. Leave everything as-is, it will be taken care of above.
                    self.inputBuffer.appendleft((a,b,i,j,k)) # return value to queue to be used after queue swap
                    self.inputFinishedFlag = True
                    self.tempI = i
                    return #since part 2 is done before part 1, this can be done.
                    
                # Merging the two queues
                if self.inputQueues[self.currentQueueNumber]: # Checks if there are any values left in the queue being merged
                    c,colnum = self.inputQueues[self.currentQueueNumber].popleft()
                    if colnum == j: # If Column number (of the value on top of the queue) and j (colnum of new data) are equal, merge and add
                        #print("merge", c + a*b, j)
                        self.inputQueues[self.helperQueueNumber].append((c + a*b,j))
                    elif colnum < j: # if colnum is less than j
                        #print("column", c, colnum)
                        self.inputQueues[self.helperQueueNumber].append((c,colnum)) # put the value in the queue into the new queue
                        self.inputBuffer.appendleft((a,b,i,j,k)) #put the value back into the input buffer, it's not getting used yet.
                    else: #if j < colnum
                        #print("add", a*b, j)
                        self.inputQueues[self.helperQueueNumber].append((a*b,j))
                        self.inputQueues[self.currentQueueNumber].appendleft((c,colnum))
                else:
                    self.inputQueues[self.helperQueueNumber].append((a*b,j)) # if the queue is empty, just add the new value.
                return
            else:
                self.partIWastedCycles += 1
            return
    
        # END SMALLMERGE
        
        
        # This is regular mode (not SmallMerge)
        if self.outputEmptyFlag and self.inputFinishedFlag:
            #print("SWAP")
            self.outputEmptyFlag = False
            self.inputFinishedFlag = False
            
            if self.endFlag2:
                self.endFlag = True
                return
            if self.endFlag1: # If we are swapping and there is no more input, we need to swap one last time so the last numbers can be cleaned up
                self.endFlag2 = True
                self.inputFinishedFlag = True
                
            #swap queue sets, set appropriate values, reset flags
            temp = self.inputQueues
            self.inputQueues = self.outputQueues
            self.outputQueues = temp
            self.helperQueueNumber = 0
            self.currentQueueNumber = 1
            self.inputQueuesLengths = [0] * self.numQueues
            
            self.prevI = self.currentI 
            self.currentI = self.tempI #remember, tempI is set when the inputFinishedFlag is set, with the new value of I
            
            # NOTE: at this point, self.inputQueues should be full of empty queues, and all of its helper variables should be reset
            # On the other hand, self.outputQueues should be full of values ready to be merged.

            
        # Part II
        if not self.outputEmptyFlag:
            
            # Check if prevI is valid!
            if self.prevI == -1:
                self.outputEmptyFlag = True
            else:
                col = -1 
                qlist = []
                # iterate through all the queues to find the indices that share the smallest column value
                for i,x in enumerate(self.outputQueues):
                    if not self.outputQueues[i]:
                        continue
                    if col == -1 or col < x[0][1]:
                        col = x[0][1]
                        qlist = [i]
                    elif col == x[0][1]:
                        qlist.append(i)
                                
                if col == -1:
                    self.outputEmptyFlag = True
                else:  
                    # now, go through those queues, popping their values and summing up their values (this is done in an adder tree)
                    # since it is done in an adder tree with a new set being added each cycle, it can be simplified to just running in one cycle
                    val = 0
                    for x in qlist:
                        val += self.outputQueues[x].popleft()[0]
                    
                    self.outputBuffer.append((val,self.prevI, col))
        else:
            self.partIIWastedCycles += 1
            

        # Part I
        if not self.inputFinishedFlag or not self.endFlag1:
            if not self.inputBuffer:
                self.inputFlag = False # If input buffer is empty, do not run
                return
            
            
            a, b, i, j, k = self.inputBuffer.popleft() #Take the input for this cycle
            
            if a == None or k != self.currentK or i != self.currentI:
                # print("FLUSHING")
                # Check that K has not been changed, if changed then flush current queue into the helper, and set it as the new helper queue
                # Remember, the "current" queue should ALWAYS be empty when swapping, because it has been flushed into the helper queue.
                # This means that we use the helper queue as the new queue to merge into, choose the shortest queue (not including helper), and merge using that and the input
                # The "shortest queue", which is being emptied, is the new "helper" for the next K change
                # Need to figure out how many cycles it will take to "flush" data from one queue to another
                
                # We also do the same operations if the new input indicates EOF. Flush queue, ensure the computation of the row has finished.
                l = len(self.inputQueues[self.currentQueueNumber])
                for x in range(l):
                    self.inputQueues[self.helperQueueNumber].append(self.inputQueues[self.currentQueueNumber].popleft())
                # adjust queue sizes
                self.inputQueuesLengths[self.helperQueueNumber] += l
                self.inputQueuesLengths[self.currentQueueNumber] = 0

                
                self.helperQueueNumber = self.currentQueueNumber
                
                # Now, we find the index of the smallest remaining queue not including the new helper queue (should be empty)
                remainder = self.inputQueuesLengths[:self.helperQueueNumber] + self.inputQueuesLengths[self.helperQueueNumber + 1:]
                self.currentQueueNumber = remainder.index(min(remainder))
                
                # if the current queue (full) is larger than the now-empty new helper queue, since it wasn't included in finding the minimum value
                # above, one must be added to the current queue number  if it is greater than the helper queue
                if self.currentQueueNumber >= self.helperQueueNumber:
                    self.currentQueueNumber += 1
                
                self.currentK = k
            
            if a == None:
                self.inputFlag = True
                self.inputFinishedFlag = True
                self.endFlag1 = True
                return # End of Processing

            if i != self.currentI:
                # If I has been changed, then set the inputFinishedFlag and change currentI, and skip. Leave everything as-is, it will be taken care of above.
                self.inputBuffer.appendleft((a,b,i,j,k)) # return value to queue to be used after queue swap
                self.inputFinishedFlag = True
                self.tempI = i
                return #since part 2 is done before part 1, this can be done.
            

                
            # The actual computation (According to paper, when k is incremented, there is still a value generated in that cycle)
            
            # Take input from input buffer
            #print("adding number")
            if self.inputQueues[self.currentQueueNumber]: # Checks if there are any values left in the queue being merged
                c,colnum = self.inputQueues[self.currentQueueNumber].popleft()
                if colnum == j: # If Column number (of the value on top of the queue) and j (colnum of new data) are equal, merge and add
                    #print("merge", c + a*b, j)
                    self.inputQueues[self.helperQueueNumber].append((c + a*b,j))
                    self.inputQueuesLengths[self.currentQueueNumber] -= 1
                elif colnum < j: # if colnum is less than j
                    #print("column", c, colnum)
                    self.inputQueues[self.helperQueueNumber].append((c,colnum)) # put the value in the queue into the new queue
                    self.inputQueuesLengths[self.currentQueueNumber] -= 1
                    self.inputBuffer.appendleft((a,b,i,j,k)) #put the value back into the input buffer, it's not getting used yet.
                else: #if j < colnum
                    #print("add", a*b, j)
                    self.inputQueues[self.helperQueueNumber].append((a*b,j))
                    self.inputQueues[self.currentQueueNumber].appendleft((c,colnum))
            else:
                self.inputQueues[self.helperQueueNumber].append((a*b,j)) # if the queue is empty, just add the new value.

            # increment length of helper queue
            self.inputQueuesLengths[self.helperQueueNumber] += 1
        else:
            self.partIWastedCycles += 1
        return

    
class SpBL:
    def __init__(self, ChannelNumber) -> None: # NEED TO FIGURE OUT SPEED OF MEMORY
        self.inputBuffer = deque([])
        self.inputFlag = False
        self.endFlag = False
        self.newAFlag = False
        self.C = ChannelNumber
        self.B_values, self.B_column_ids, self.B_row_lengths, self.B_row_pointers = (None,None,None,None) # Stored in C2SR format
        self.PE = None
        self.currentValueBuffer = deque([]) # This is a buffer containing all of the stuff to be sent to the PE one by one.
        self.currentColumnBuffer = deque([])
        self.a = 0
        self.i = 0
        self.k = 0
        
        self.memory = None
        self.memoryWastedCycles = 0
        self.MemoryUsage = 0
        self.freeBytes = 0
        
    def __str__(self) -> str:
        return "i: " + str(self.i) + ", k:" + str(self.k) + ", ValueBuffer: " + str(self.currentValueBuffer) + ", ColumnBuffer: " + str(self.currentColumnBuffer)
        
    def setNext(self,p: "PE"):
        self.PE = p
    
    def setMemory(self,m: "Memory"):
        self.memory = m

    def input(self, a, i, k):
        self.inputBuffer.append((a,i,k))
        self.inputFlag = True
        
    def loadMatrixB(self, BMatrix):
        self.B_values, self.B_column_ids, self.B_row_lengths, self.B_row_pointers = BMatrix
        
    
    def cycle(self):
        if not self.inputFlag or self.endFlag:
            return
        
        if not self.currentValueBuffer:
            
            if self.inputBuffer:
                self.a, self.i, self.k = self.inputBuffer.popleft()
            else:
                self.inputFlag = False
                return
            
            if self.a == None:
                self.endFlag = True
                self.PE.input(None,None,None,None,None)
                return
        
            # Because we are loading from matrix B, We need to check where each row is (channel and location) using row lengths and pointers
            self.currentValueBuffer = deque(self.B_values[self.B_row_pointers[self.k][0]][self.B_row_pointers[self.k][1]:self.B_row_pointers[self.k][1] + self.B_row_lengths[self.k][1]])
            self.currentColumnBuffer = deque(self.B_column_ids[self.B_row_pointers[self.k][0]][self.B_row_pointers[self.k][1]:self.B_row_pointers[self.k][1] + self.B_row_lengths[self.k][1]])
            
            numBytes = len(self.currentColumnBuffer) * 4 + len(self.currentValueBuffer) * 4
            if numBytes == 0:
                return
            # 64 bytes are requested in vector fashion (channel interleaving size = 64) - since we use C2SR, no data is wasted.
            for x in range(0,numBytes//64+1):
                ms = min(64,numBytes)
                self.memory.requestData(self,self.B_row_pointers[self.k][0],ms)
                self.MemoryUsage += ms
        else:
            if self.freeBytes < 8:
                self.memoryWastedCycles += 1
                return
            # 8 bytes for one set of column, value
            self.freeBytes -= 8
            # in the else statement because if a row is empty we could pop from an empty queue (also fetching data should take at least one cycle to queue it lol)
            b = self.currentValueBuffer.popleft()
            j = self.currentColumnBuffer.popleft()
            self.PE.input(self.a,b,self.i,j,self.k)        
        


    
    
class SpAL:
    def __init__(self, ChannelNumber, numChannels) -> None:
        self.endFlag = False
        self.C = ChannelNumber
        self.numChannels = numChannels
        self.A_values, self.A_column_ids, self.A_row_lengths, self.A_row_pointers = (None, None, None, None)
        self.SpBL = None
        self.memory = None
        self.i = ChannelNumber # start at 0, and add numChannels every time
        self.currentValueBuffer = None
        self.currentColumnBuffer = None
        self.memoryWastedCycles = 0
        self.MemoryUsage = 0
        self.freeBytes = 0
        
    def setNext(self,b: "SpBL"):
        self.SpBL = b
    
    def setMemory(self,m: "Memory"):
        self.memory = m
    
    def loadMatrixA(self, AMatrix):
        self.A_values, self.A_column_ids, self.A_row_lengths, self.A_row_pointers = AMatrix # Stored in C2SR format
        self.currentValueBuffer = deque(self.A_values[self.C][self.A_row_pointers[self.i][1]:self.A_row_pointers[self.i][1] + self.A_row_lengths[self.i][1]]) # First set of values (one row of A)
        self.currentColumnBuffer = deque(self.A_column_ids[self.C][self.A_row_pointers[self.i][1]:self.A_row_pointers[self.i][1] + self.A_row_lengths[self.i][1]]) # First set of column numbes (one row of A)
        numBytes = len(self.currentColumnBuffer)*4 + len(self.currentValueBuffer)*4 # 32-bit values and column numbers
        self.memory.requestData(self,self.C,numBytes)
        self.MemoryUsage += numBytes
    
    def __str__(self) -> str:
        return "i: " + str(self.i) + ", ValueBuffer: " + str(self.currentValueBuffer) + ", ColumnBuffer: " + str(self.currentColumnBuffer)
    
    def cycle(self):

        if self.endFlag:
            return
        
        if not self.currentValueBuffer:
            self.i += self.numChannels
            if self.i >= len(self.A_row_pointers):
                self.SpBL.input(None,None,None)
                self.endFlag = True
                return
            self.currentValueBuffer = deque(self.A_values[self.C][self.A_row_pointers[self.i][1]:self.A_row_pointers[self.i][1] + self.A_row_lengths[self.i][1]])
            self.currentColumnBuffer = deque(self.A_column_ids[self.C][self.A_row_pointers[self.i][1]:self.A_row_pointers[self.i][1] + self.A_row_lengths[self.i][1]])
            numBytes = len(self.currentColumnBuffer)*4 + len(self.currentValueBuffer)*4 # 32-bit values and column numbers
            if numBytes == 0:
                return
            for x in range(0,numBytes//64+1):
                ms = min(64,numBytes)
                self.memory.requestData(self,self.A_row_pointers[self.i][0],ms)
                self.MemoryUsage += ms
                numBytes -= 64
        else:
            if self.freeBytes < 8:
                self.memoryWastedCycles += 1
                return
            # 8 bytes for one set of column, value
            self.freeBytes -= 8
            k = self.currentColumnBuffer.popleft()
            a = self.currentValueBuffer.popleft() 
            self.SpBL.input(a,self.i,k)
            
class Memory:
    def __init__(self, numChannels, peakBandwidthPerChannel) -> None:
        self.numChannels = numChannels
        self.channelQueues = [deque() for _ in range(numChannels)]
        self.channelCurrent = [(None,0,0) for _ in range(numChannels)]
        self.BPC = peakBandwidthPerChannel
    
    def requestData(self,loader, channel, amount):
        # numCycles includes the row and column latencies
        self.channelQueues[channel].append([loader,amount]) # (loader, number of bytes, memory latency)
    
    def cycle(self):
        for i,channel in enumerate(self.channelCurrent):
            if channel[1] > 0:
                sent = min(channel[1],self.BPC)
                channel[0].freeBytes += sent # "send" the bytes over to the PE.
                channel[1] -= sent
            elif channel[1] <= 0 and self.channelQueues[i]:
                self.channelCurrent[i] = self.channelQueues[i].popleft()

