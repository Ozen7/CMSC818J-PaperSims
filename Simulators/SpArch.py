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
    
    