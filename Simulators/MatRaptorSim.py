import numpy as np
# Step 1: Convert CSR to C2SR - C2SR basically just makes sure that all the data a PE reads at a given time is useful.

# Step 2: Develop functionality of a single PE assigned to each row.


class Controller:
    def __init__(self,numPEs,queuesize,M1,M2):
        pass

class PE:
    def __init__(self,queuesize):
        self.queuelist = np.zeros((9,queuesize,2))
        self.helperqueue = np.zeros(queuesize)
        self.queuelengths = np.zeros(10)        
        
    def merge(self,val,queuenum):
        pass
        
    # the paper was confusing on whether it included K or not in its simulation.
    def compute(self, a, b,i,j,k):
        if self.kprev == None or self.kprev != k:
            self.kprev = k
            #This is not part of the original design, but its the only way I could figure out how to do it
            self.currentqueue = np.argmin(self.queuelengths)
                
        
        c = a * b
        self.merge(c,self.currentqueue)
        
        
        
            
    
    
    
    
    
    

    
    

# Step 3: add multiple Processing elements, and optimize.

# Step 4: add cycle counters as well as memory bandwidth counters



