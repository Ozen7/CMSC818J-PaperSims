import numpy as np
from scipy.sparse import csr_matrix



#Matrix Condensing
class Afetcher:
    def __init__(self, A) -> None:
        self.A = csr_matrix(A)
        data = self.A.data
        indices = self.A.indices
        self.originalcols = np.copy(A.indices)
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
    
    #returns a set of values that are in the same row of the current compressed column.
    def nextRowSet(self):
        condenseSize = 3
        if self.rowPointer > self.A.shape[0]:
            self.rowPointer = 1
            self.colPointer -= condenseSize
        minval = 0 if self.colPointer - condenseSize <= 0 else self.colPointer-condenseSize
        rp0 = indptr[self.rowPointer-1]
        rp1 = indptr[self.rowPointer]
        mask = np.in1d(self.A.indices[rp0: rp1],np.arange(minval, self.colPointer))
        while not mask.any():
            if self.rowPointer >= self.A.shape[0]:
                self.rowPointer = 1
                self.colPointer -= condenseSize
                if self.colPointer <= 0:
                    return None
                minval = 0 if self.colPointer - condenseSize <= 0 else self.colPointer-condenseSize
            else:
                self.rowPointer += 1
            rp0 = indptr[self.rowPointer-1]
            rp1 = indptr[self.rowPointer]
            mask = np.in1d(self.A.indices[rp0: rp1],np.arange(minval, self.colPointer))
        #print(mask)
        #print(self.originalcols)
        #print(rp0)
        #print(rp1)
        #print(self.A.indices)
        #print(self.originalcols)
        # [original(and compressed) row, which columns within the compressed array, original columns, data]
        retval = [self.rowPointer-1, self.A.indices[rp0:rp1][mask],  self.originalcols[rp0:rp1][mask],self.A.data[rp0:rp1][mask]]
        
        self.rowPointer += 1
        return retval

# Row Prefetcher
class Bprefetcher:
    def __init__(self) -> None:
        pass


class MultiplyAndMerge:
    def __init__(self, A, B) -> None:
        pass
    
        
        

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
print(f.nextRowSet())


    

#Huffman Tree Scheduler, Comparator Array Based Merger.
# NOTE: The pipelined multiply and merge can actually just be done manually 
# since it only takes one clock cycle

# Partial Marices are represented in COO format Sorted by row index then col index
# Instead of just going through and merging the values normally,
class mergeUnit:
    def __init__(self) -> None:
        pass
    
def outerProduct(M1, M2):
    pass

