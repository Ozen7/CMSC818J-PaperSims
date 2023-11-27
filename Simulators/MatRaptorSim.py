import numpy as np
from collections import deque
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

import scipy.io as sio

sparseMat = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')


# Step 1: Convert CSR to C2SR - C2SR basically just makes sure that all the data a PE reads at a given time is useful.
def csr_to_c2sr(data, indices, indptr, num_channels):
    # Create arrays to store C2SR format data
    c2sr_values = [[] for x in range(num_channels)]
    c2sr_column_ids = [[] for x in range(num_channels)]
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
"""
# Step 2: Develop functionality of a single PE assigned to each row.

class PE:
    def __init__(self,queuesize,queuenum):
        self.queuenum = queuenum
        self.queuesize = queuesize
        self.queuelist = [deque(maxlen=queuesize) for x in range(0,queuenum)]
        self.helperqueue = 0
        self.nomerge = True
        self.kprev = -1
    
    def reset(self):
        self.queuelist = [deque(maxlen=self.queuesize) for x in range(0,self.queuenum)]
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
                for x in range(0,len(self.queuelist[self.currentqueue])):
                    self.queuelist[self.helperqueue].append(self.queuelist[self.currentqueue].popleft())
                    
                assert(len(self.queuelist[self.currentqueue]) == 0)
                self.helperqueue = self.currentqueue
                
            self.kprev = k
            
            #This is not part of the original design, but its the only way I could figure out how to do it
            adjustedqueues = self.queuelist[:self.helperqueue] + self.queuelist[self.helperqueue + 1:]

            self.currentqueue = self.queuelist.index(min(adjustedqueues, key=len, default=None))
            
            if len(self.queuelist[self.currentqueue]) == 0:
                self.nomerge = True
            else:
                self.nomerge = False
                
        
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
            
            mins = np.where(firstcols == firstcols.min())[0]
            indices.append(self.queuelist[mins[0]][0][1])

            for num in mins:
                total += self.queuelist[num].popleft()[0]
            data.append(total)
        self.reset()
        return (data, indices)
        
            
                
        
        
# Step 3: Develop a way to load information in parallel to each PE, splitting based on row. 
class Controller:
    def __init__(self,A,B,numchannels,queuelengths,queuenum):
        #Both of these are converted to C2SR format.
        self.A_values, self.A_column_ids, self.A_row_lengths, self.A_row_pointers = csr_to_c2sr(A.data, A.indices,A.indptr,numchannels)
        self.B_values, self.B_column_ids, self.B_row_lengths, self.B_row_pointers = csr_to_c2sr(B.data, B.indices,B.indptr,numchannels)
        
        self.outputshape = (A.shape[0], B.shape[1])
        
        self.pes = [PE(queuesize=queuelengths,queuenum=queuenum) for x in range(numchannels)]
            
        self.numchannels = numchannels
    
    def compute(self):
        indices = []
        indptr = []
        data = []
        for row in range(len(self.A_row_lengths)):
            channel = row % self.numchannels
            PE = self.pes[channel]
            rowstart = self.A_row_pointers[row][1]
            
            for a in range(self.A_row_lengths[row][1]):
                adata = self.A_values[channel][rowstart + a]
                acol = self.A_column_ids[channel][rowstart + a]
                (browchannel, browstart) = self.B_row_pointers[acol]
                
                for b in range(self.B_row_lengths[acol][1]):
                    bdata = self.B_values[browchannel][browstart + b]
                    bcol = self.B_column_ids[browchannel][browstart + b]
                    PE.compute(adata,bdata,row,acol,bcol)
            (d, i) = PE.accumulateandoutput()
            indptr.append(len(data))
            indices += i
            data += d
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

#controller = Controller(csr_matrix(i1),csr_matrix(i2),10,10,100)
controller = Controller(csr_matrix(sparseMat),csr_matrix(sparseMat),10,100000,100)

out = controller.compute()


print(out.toarray())
print(np.matmul(sparseMat.toarray(),sparseMat.toarray()))
#print(np.matmul(i1.toarray(),i2.toarray()))
print(np.equal(out.toarray(), np.matmul(sparseMat.toarray(),sparseMat.toarray())))
print(np.allclose(out.toarray(),np.matmul(sparseMat.toarray(),sparseMat.toarray()),rtol=0.000001))

# pe = PE(50)

# A
# [1, 0, 2, 3]
# 
# B
# [1, 0, 0, 1]
# [0, 0, 0, 0]
# [1, 0, 1, 0]
# [1, 1, 0, 1]
#
# C (should equal):
# [6, 3, 2, 4] 
#pe.compute(1,1,0,0,0)
#pe.compute(1,1,0,0,3)

#pe.compute(2,1,0,2,0)
#pe.compute(2,1,0,2,2)

#pe.compute(3,1,0,3,0)
#pe.compute(3,1,0,3,1)
#pe.compute(3,1,0,3,3)

#print(pe.queuelist)


# Step 4: Add true parallelism to the code

# Step 5: add cycle counters as well as memory bandwidth counters



