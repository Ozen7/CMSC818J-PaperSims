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
    
    def __init__(self, numQueues, PENumber) -> None:
        self.inputBuffer = deque([])
        self.outputBuffer = deque([])
        
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
            
            if a == None or k != self.currentK:
                #print("FLUSHING")
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
        
    def __str__(self) -> str:
        return "i: " + str(self.i) + ", k:" + str(self.k) + ", ValueBuffer: " + str(self.currentValueBuffer) + ", ColumnBuffer: " + str(self.currentColumnBuffer)
        
    def setNext(self,p: "PE"):
        self.PE = p
    
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
        else:
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
        self.i = ChannelNumber # start at 0, and add numChannels every time
        self.currentValueBuffer = None
        self.currentColumnBuffer = None
        
    def setNext(self,b: "SpBL"):
        self.SpBL = b
    
    def loadMatrixA(self, AMatrix):
        self.A_values, self.A_column_ids, self.A_row_lengths, self.A_row_pointers = AMatrix # Stored in C2SR format
        self.currentValueBuffer = deque(self.A_values[self.C][self.A_row_pointers[self.i][1]:self.A_row_pointers[self.i][1] + self.A_row_lengths[self.i][1]]) # First set of values (one row of A)
        self.currentColumnBuffer = deque(self.A_column_ids[self.C][self.A_row_pointers[self.i][1]:self.A_row_pointers[self.i][1] + self.A_row_lengths[self.i][1]]) # First set of column numbes (one row of A)
    
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

        else:
            # this is in the else because a row could be empty!
            #print("LEN OF COLUMN", self.C, len(self.currentColumnBuffer))
            #print("LEN OF VALUE", self.C, len(self.currentValueBuffer))
            k = self.currentColumnBuffer.popleft()
            a = self.currentValueBuffer.popleft() 
            self.SpBL.input(a,self.i,k)

"""
# Example usage
# data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# indices = [0, 2, 4, 1, 3, 4, 0, 2, 3]
# indptr = [0, 2, 4, 7, 9]
datacsr = csr_matrix(sparsematrix)
data = datacsr.data
indices = datacsr.indices
indptr = datacsr.indptr

num_channels = 3

c2sr_values, c2sr_column_ids, c2sr_row_length, c2sr_row_pointer = csr_to_c2sr(data, indices, indptr, num_channels)

# Print the results
assert(len(c2sr_values) == len(c2sr_column_ids))

print("C2SR Values:", c2sr_values)
print("C2SR Column IDs:", c2sr_column_ids)
print("C2SR Row Lengths:", c2sr_row_length)
print("C2SR Row Pointers:", c2sr_row_pointer)


# Step 2: Develop functionality of a single PE assigned to each row.

class PE:
    def __init__(self,queuesize,queuenum):
        self.phase1Cycles = 0
        self.phase2Cycles = 0
        self.queuenum = queuenum
        self.queuesize = queuesize
        self.queuelist = [deque(maxlen=queuesize) for _ in range(0,queuenum)]
        self.helperqueue = 0
        self.nomerge = True
        self.kprev = -1
    
    def getNumCycles(self):
        return (self.phase1Cycles, self.phase2Cycles)
    
    def reset(self):
        self.queuelist = [deque(maxlen=self.queuesize) for _ in range(0,self.queuenum)]
        self.phase1Cycles = 0
        self.phase2Cycles = 0
        self.helperqueue = 0
        self.nomerge = True
        self.kprev = -1
        
    def merge(self,val,i,j):
        if self.nomerge:
            self.queuelist[self.currentqueue].append((val,j))
        else:
            while len(self.queuelist[self.currentqueue]) != 0 and self.queuelist[self.currentqueue][0][1] < j:
                # loop through queue until we find a col value greater than or equal to that of our new input
                self.queuelist[self.helperqueue].append(self.queuelist[self.currentqueue].popleft())
            
            #If they're equal in column, add them up instead.
            if len(self.queuelist[self.currentqueue]) != 0 and self.queuelist[self.currentqueue][0][1] == j:
                self.queuelist[self.helperqueue].append((val + self.queuelist[self.currentqueue].popleft()[0],j))
            else:
                self.queuelist[self.helperqueue].append((val,j))
                
        
    # the paper was confusing on whether it included K or not in its simulation.
    def compute(self, a, b, i, k, j):
        
        # if there is no prior k value, or the prior k value is not the same as the current one 
        # (we are operating on a new row in B)
        if self.kprev != k:
            # if we DID merge in the previous vector with another queue, we now flush the remaining values into the helper queue.
            # the now-empty queue we merged into is now the helper queue
            if not self.nomerge:
                for _ in range(0,len(self.queuelist[self.currentqueue])):
                    self.queuelist[self.helperqueue].append(self.queuelist[self.currentqueue].popleft())

                
                assert(len(self.queuelist[self.currentqueue]) == 0)
                self.helperqueue = self.currentqueue
                
            self.kprev = k
            
            #This is not part of the original design, but its the only way I could figure out how to do it
            adjustedqueues = self.queuelist[:self.helperqueue] + self.queuelist[self.helperqueue + 1:]
            newqueueidx = adjustedqueues.index(min(adjustedqueues, key=len, default=None))
            if newqueueidx >= self.helperqueue:
                newqueueidx += 1
            self.currentqueue = newqueueidx
            
            if len(self.queuelist[self.currentqueue]) == 0:
                self.nomerge = True
            else:
                self.nomerge = False
                
        self.phase1Cycles += 1 # one cycle for each multiplication, the merging and accumulation is masked by another multiplication
        c = a * b
        self.merge(c,i,j)
    
    def accumulateandoutput(self):
        data = []
        indices = []
        while True:
            firstcols = np.zeros(self.queuenum)
            for x in range(self.queuenum):
                if self.queuelist[x]:
                    firstcols[x] = self.queuelist[x][0][1]
                else:
                    firstcols[x] =  2**30
            total = 0
        
            for x in firstcols:
                if x != 2**30:
                    break
            else:
                break
            
            # find minimum values across all queues and accumulate them. according to the paper, this ONLY takes one cycle. 
            # "In the case when multiple queues have the same minimum column index at the top of the queue, all such queues are popped and
            # the sum of popped elements is streamed to the main memory"
            # Note that there is a dedicated portion that keeps track of the minimum column IDs as well as an adder tree, so this isn't too 
            # far fetched, even if the addition should probably technically take more than one cycle.
            self.phase2Cycles += 1
            mins = np.where(firstcols == firstcols.min())[0]
            indices.append(self.queuelist[mins[0]][0][1])

            for num in mins:
                total += self.queuelist[num].popleft()[0]
            data.append(total)
        return (data, indices)
        
            
                
        
        
# Step 3: Develop a way to load information in parallel to each PE, splitting based on row. 
class Controller:
    def __init__(self,A,B,numchannels,queuelengths,queuenum):
        #Both of these are converted to C2SR format.
        self.A_values, self.A_column_ids, self.A_row_lengths, self.A_row_pointers = csr_to_c2sr(A.data, A.indices,A.indptr,numchannels)
        self.B_values, self.B_column_ids, self.B_row_lengths, self.B_row_pointers = csr_to_c2sr(B.data, B.indices,B.indptr,numchannels)

        self.outputshape = (A.shape[0], B.shape[1])
        
        self.pes = [PE(queuesize=queuelengths,queuenum=queuenum) for _ in range(numchannels)]
            
        self.numchannels = numchannels
        
        self.peNumCycles = [[0,0] for _ in range(self.numchannels)]
    
    def obtainMaxCycles(self):
        return max(map(lambda x : x[0] + x[1], self.peNumCycles))
    
    def compute(self):
        global total_memory_usage
        indices = []
        indptr = []
        data = []
        channel = 0
        finaladd = 0
        for row in range(len(self.A_row_lengths)):
            total_memory_usage += 2 # read from A: (length, row pointer)
            channel = row % self.numchannels
            PE = self.pes[channel]
            rowstart = self.A_row_pointers[row][1]
            
            for a in range(self.A_row_lengths[row][1]):
                total_memory_usage += 2 # read from A: (value, column id)
                adata = self.A_values[channel][rowstart + a]
                acol = self.A_column_ids[channel][rowstart + a]
                (browchannel, browstart) = self.B_row_pointers[acol]
                total_memory_usage += 2 # read from B: (length, row pointer)
                
                for b in range(self.B_row_lengths[acol][1]):
                    total_memory_usage += 2 # read from B: (value, column id)
                    bdata = self.B_values[browchannel][browstart + b]
                    bcol = self.B_column_ids[browchannel][browstart + b]
                    PE.compute(adata,bdata,row,acol,bcol)
            (d, i) = PE.accumulateandoutput()
            (phase1,phase2) = PE.getNumCycles()
            self.peNumCycles[channel][0] += max(phase1,self.peNumCycles[channel][1])
            self.peNumCycles[channel][1] = phase2
            finaladd = phase2
            PE.reset()
            indptr.append(len(data))
            indices += i
            data += d
        self.peNumCycles[channel][0] += finaladd
        indptr.append(len(data))
        
        if indptr != []:
            return csr_matrix((data,indices,indptr), shape=(self.outputshape))
        else:
            return 0
            

#A = [1,0,2,3]
#B = [[1,0,0,1],[0,0,0,0],[1,0,1,0],[1,1,0,1]]
#controller = Controller(csr_matrix(A),csr_matrix(B),2,50,3)

gen = np.random.default_rng()
data1 = gen.integers(0,10,5)
row1 = gen.integers(0,5,5)
col1 = gen.integers(0,5,5)

data2 = gen.integers(0,10,5)
row2 = gen.integers(0,5,5)
col2 = gen.integers(0,5,5)
i1 = coo_matrix((data1, (row1, col1)), shape=(5, 5))
i2 = coo_matrix((data2, (row2, col2)), shape=(5, 5))

'''
sparseMat = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')
controller = Controller(csr_matrix(sparseMat),csr_matrix(sparseMat),10,100000,100)

out = controller.compute()


print(out.toarray())
print(np.matmul(sparseMat.toarray(),sparseMat.toarray()))
#print(np.matmul(i1.toarray(),i2.toarray()))
print(np.equal(out.toarray(), np.matmul(sparseMat.toarray(),sparseMat.toarray())))
print(np.allclose(out.toarray(),np.matmul(sparseMat.toarray(),sparseMat.toarray()),rtol=0.000001))

'''
#controller = Controller(csr_matrix(i1),csr_matrix(i2),10,10,100)

# input matrices are in CSR format
def run_matraptor(matrix1,matrix2):
    csr_m1 = csr_matrix(matrix1)
    csr_m2 = csr_matrix(matrix2)

    controller = Controller(csr_m1,csr_m2,10,100000,100)

    out = controller.compute()
        
    global total_memory_usage

    print("number of cycles:", controller.obtainMaxCycles())
    print("Data Pulled from DRAM:", total_memory_usage)
    if matrix1.size >= 10000 or matrix2.size >= 10000:
        print("matrices too large to verify")
    #else:
    #    trueval = matrix1 @ matrix2
    #    print("Verify that sparse multiplication is correct: ", np.allclose(out.toarray(),trueval,rtol=0.000001))
    total_memory_usage = 0
"""
