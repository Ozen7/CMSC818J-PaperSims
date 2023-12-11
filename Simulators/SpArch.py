import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from collections import deque


#Matrix Condensing
class Afetcher:
    def __init__(self, A) -> None:
        self.condenseSize = 63
        self.A = csr_matrix(A)
        data = self.A.data
        indices = self.A.indices
        self.originalcols = np.copy(self.A.indices)
        indptr = self.A.indptr
        
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
        #print(mask)
        #print(self.originalcols)
        #print(rp0)
        #print(rp1)
        #print(self.A.indices)
        #print(self.originalcols)
        # [original(and compressed) row, which columns within the compressed array, original columns, data]
        retval = [self.rowPointer-1, compressedcols,  self.originalcols[rp0:rp1][mask],self.A.data[rp0:rp1][mask]]
        
        self.rowPointer += 1
        return retval

# Row Prefetcher
class Bprefetcher:
    def __init__(self,B) -> None:
        self.B = csr_matrix(B)
    def fetch(self, row):
        begin = self.B.indptr[row]
        end = self.B.indptr[row+1]
        return [self.B.data[begin:end],self.B.indices[begin:end]]


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
        
            

global truepartial
truepartial = []
class MultiplyAndMerge:
    def __init__(self, A, B) -> None:
        global truepartial
        self.aFetcher = Afetcher(A)
        self.bPrefetcher = Bprefetcher(B)
        self.partials = [deque() for _ in range(64)]
        for x in truepartial:
            self.partials[63].append(x)
        self.endflag = False
        
    def MultiplyColumn(self):
        
        rs = self.aFetcher.nextRowSet()
        while rs != None:
            for x in range(len(rs[3])):
                brow = np.array(self.bPrefetcher.fetch(rs[2][x]))
                datapoints = brow[0] * rs[3][x]
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


    
    
    
        
        
'''
data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
indices = [0, 2, 4, 1, 3, 4, 0, 2, 3]
indptr = [0, 2, 4, 7, 9]
datacsr = csr_matrix((data,indices,indptr))
datacsr.shape
print(datacsr.toarray())
f = Afetcher(datacsr)
print(f)
print(f.nextRowSet())
print(f.nextRowSet())
print(f.nextRowSet())
print(f.nextRowSet())
print(f.nextRowSet())
f.resetRows()
print(f.nextRowSet())
print(f.nextRowSet())
print(f.nextRowSet())
print(f.nextRowSet())
print(f.nextRowSet())


'''
gen = np.random.default_rng()
data1 = gen.integers(1,10,10000)
row1 = gen.integers(0,1000,10000)
col1 = gen.integers(0,1000,10000)

data2 = gen.integers(1,10,10000)
row2 = gen.integers(0,1000,10000)
col2 = gen.integers(0,1000,10000)
i1 = coo_matrix((data1, (row1, col1)), shape=(1000, 1000))
i2 = coo_matrix((data2, (row2, col2)), shape=(1000, 1000))

print(i1.toarray())
print(i2.toarray())

merger = MultiplyAndMerge(i1,i2)
merger.MultiplyColumn()
merger.MergePartials()

data = []
i = []
j = []

for x in truepartial:
    i.append(x[0])
    j.append(x[1])
    data.append(x[2])
    
true = np.matmul(i1.toarray(),i2.toarray())
m = coo_matrix((data, (i, j)), shape= true.shape)
print(m.toarray())
print(true)

print(np.allclose(m.toarray(),true,0.0001))