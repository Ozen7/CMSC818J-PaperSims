import numpy as np
from collections import deque
# Step 1: Convert CSR to C2SR - C2SR basically just makes sure that all the data a PE reads at a given time is useful.

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
        
        
pe = PE(50)

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


            
    
    
    
    
    
    

    
    

# Step 3: add multiple Processing elements, and optimize.

# Step 4: add cycle counters as well as memory bandwidth counters



