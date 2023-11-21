from typing import Callable
import numpy as np
import scipy.io as sio

data = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')

print(type(data))
class Scanner:

    def __init__(self) -> None:
        self.stack = []
    
    def iterate(self) -> int:
        return self.stack.pop(0)

    def readNew(self, input) -> None:
        self.flush()
        self.stack = input
        
    def isEmpty(self) -> bool:
        return self.stack == []

    def flush(self) -> None:
        self.stack = []
        
    def peek(self) -> int:
        return self.stack[0]

    

class Intersect:

        
    def reset(self) -> None:
        self.A = None
        self.B = None
        self.op = None
        self.output = np.array([])
    
    def set(self,scanner1: Scanner, scanner2: Scanner, operation: Callable[[int,int],int] ) -> None:
        self.A = scanner1
        self.B = scanner2
        self.op = operation
        self.output = np.array([])
    
    
    def getScanner1(self) -> Scanner:
        return self.A
    
    def getScanner2(self) -> Scanner:
        return self.B

    def getOutput(self):
        return self.output
    
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
                self.output = np.append(self.output,self.A.iterate())
                
intersector = Intersect()
scanner1 = Scanner()
scanner2 = Scanner()

    
scanner1.readNew(sorted(list(np.random.randint(0,10000,800))))
scanner2.readNew(sorted(list(np.random.randint(0,10000,800))))

intersector.set(scanner1,scanner2,lambda x, y: x*y)

intersector.Intersect()

print(intersector.output)
    