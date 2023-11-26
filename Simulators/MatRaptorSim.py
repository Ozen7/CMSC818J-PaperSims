import numpy as np
from collections import deque
# Step 1: Convert CSR to C2SR - C2SR basically just makes sure that all the data a PE reads at a given time is useful.
def csr_to_c2sr(data, indices, indptr, num_channels):
    # Create arrays to store C2SR format data
    c2sr_values = [[] for x in range(num_channels)]
    c2sr_column_ids = [[] for x in range(num_channels)]
    c2sr_row_length = [[] for x in range(num_channels)]
    c2sr_row_pointer = [[] for x in range(num_channels)]

    # Loop through each row in CSR format
    for i in range(len(indptr) - 1):
        # Get the start and end indices for the current row in CSR format
        start_idx = indptr[i]
        end_idx = indptr[i + 1]

        # Get the channel assigned to the current row
        channel = i % num_channels
        
        # Append the row length and row pointer for the current row in C2SR format
        c2sr_row_length[channel].append(end_idx - start_idx)
        c2sr_row_pointer[channel].append(len(c2sr_values[channel]))

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


class Controller:
    def __init__(self,numPEs,queuesize,M1,M2):
        pass

class PE:
    def __init__(self,queuesize):
        self.queuelist = [deque(maxlen=queuesize) for x in range(0,3)]
        self.helperqueue = 2
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
        print(self.queuelist)

                
        
    # the paper was confusing on whether it included K or not in its simulation.
    def compute(self, a, b, i, k, j):
        
        # if there is no prior k value, or the prior k value is not the same as the current one 
        # (we are operating on a new row in B)
        if self.kprev != k:
            # if we DID merge in the previous vector with another queue, we now flush the remaining values into the helper queue.
            # the now-empty queue we merged into is now the helper queue
            if not self.nomerge:
                for x in range(0,len(self.queuelist[self.currentqueue])):
                    self.helperqueue.append(self.queuelist[self.currentqueue].popleft())
                    
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
        
# Step 3: Develop a way to load information in parallel to each PE, splitting based on row. 
class Loader:
    
    def __init__(self,A,B):
        
        

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
pe.compute(1,1,0,0,0)
pe.compute(1,1,0,0,3)

pe.compute(2,1,0,2,0)
pe.compute(2,1,0,2,2)

pe.compute(3,1,0,3,0)
pe.compute(3,1,0,3,1)
pe.compute(3,1,0,3,3)

print(pe.queuelist)


# Step 4: Add true parallelism to the code

# Step 5: add cycle counters as well as memory bandwidth counters



