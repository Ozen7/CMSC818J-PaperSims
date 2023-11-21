from typing import Callable
import numpy as np
import scipy.io as sio
from scipy.sparse import csc_array
from scipy.sparse import csr_array


data = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')
cscData = csc_array(data)
csrData = csr_array(data)


# Important note: This simulator ONLY works for 2-tensors, as that is what I need for my own testing purposes.
# As such, it makes use of a COO compression format as opposed to CSF or some other higher-dimensional compression format.
print(type(data))
print(cscData.indices)


        
class Scanner:
    
    metadata = []
    level = 0
    parent_node = 0 
    
    # Metadata is simply an array of coordinates. By repeatedly intersecting two coordinate streams, we can obtain the 
    # final pairing of operations necessary to get the output.
    def loadData(self, metadata, level, parent_node) -> None:
        self.metadata = metadata # to make pop() faster
        self.metadata.reverse()
        self.level = level
        self.parent_node = parent_node
    
    def iterate(self) -> int:
        return self.metadata.pop()

    def isEmpty(self) -> bool:
        return self.metadata == []

    def flush(self) -> None:
        self.metadata = []
        
    def peek(self) -> int:
        return self.metadata[-1]
    
    def getLevel(self) -> int:
        return self.level
    
    def getParentNode(self) -> int:
        return self.parent_node
        

    

class Intersect:  
    def reset(self) -> None:
        self.A = None
        self.B = None
        self.op = None
        self.output = np.array([])
    
    def set(self,scanner1: Scanner, scanner2: Scanner ) -> None:
        self.A = scanner1
        self.B = scanner2
        self.output = np.array([])
    
    def getScanner1(self) -> Scanner:
        return self.A
    
    def getScanner2(self) -> Scanner:
        return self.B

    def getOutput(self):
        #the levels of both should be the same, since its output stationary.
        return (self.output,self.A.getParentNode(),self.B.getParentNode(), self.A.getLevel())
    
    def Intersect(self):
        while True:
            if self.A.isEmpty() or self.B.isEmpty():
                self.A.flush()
                self.B.flush()
                break
            elif self.A.peek() < self.B.peek():
                self.A.iterate();
            elif self.A.peek() > self.B.peek():
                self.B.iterate();
            else:
                self.B.iterate()
                self.output = np.append(self.output,(self.A.iterate()))
                
                
# this all will be part of the controller class.
intersector = Intersect()
scanner1 = Scanner()
scanner2 = Scanner()

gen = np.random.default_rng()
input1 =gen.integers(0,10000,800)
input2 = gen.integers(0,10000,800)
input1.sort()
input2.sort()

scanner1.loadData(list(input1),0,0)
scanner2.loadData(list(input2),0,0)

intersector.set(scanner1,scanner2)

intersector.Intersect()

for x in intersector.getOutput():
    
    #generate new random data for each "column" and "row"
    input1 =gen.integers(0,10000,800)
    input2 = gen.integers(0,10000,800)
    input1.sort()
    input2.sort()

    intersector.reset()
    scanner1.loadData()
    intersector.set()

print(intersector.getOutput())
    