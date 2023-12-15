from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
import scipy
import numpy as np
import math

        

class PE:
    def __init__(self,row,col,systolicArray) -> None:
        self.row = row
        self.col = col
        self.systolic = systolicArray
        self.value = 0
        self.column = -1
        if self.col == 0:
            self.Ain = 0
        else:
            self.Ain = None
        self.xval = 0
        self.xidx = 0 
        self.xend = False
        self.val = 0
        self.endflag = False
    
    def reset(self):
        self.value = 0
        self.column = -1
        if self.col == 0:
            self.Ain = 0
        else:
            self.Ain = None
        self.xval = 0
        self.xidx = 0 
        self.xend = False
        self.val = 0
        self.endflag = False
        

    def __str__(self) -> str:
        return str(self.value)

    def loadval(self,value,column):
        self.value = value
        self.column = column
    
    def inputLeft(self,Ain):
        self.Ain = Ain
    
    def inputTop(self,xval,xidx,xend):
        self.xval = xval
        self.xidx = xidx
        self.xend = xend
    
    def accumulateAndOutputRight(self):
        if self.endflag:
            return
        self.endflag = True
        acc = self.Ain + self.val
        if self.col + 1 >= 128:
            self.systolic.accumulate[self.row] += acc
        else:
            self.systolic.PEarray[self.row][self.col+1].inputLeft(acc)
            
    def cycle(self):
        if self.xend: #accumulate mode - wait for Ain
            if self.Ain != None: 
                self.accumulateAndOutputRight()
        elif self.xidx == self.column and self.value != 0 and self.xval != 0: # Latch Mode - save input
            self.val = self.xval * self.value
            
        # hold mode / bypass mode
        if self.row < 127:
            self.systolic.PEarray[self.row+1][self.col].inputTop(self.xval,self.xidx,self.xend)

        return self.endflag
        
class SystolicArray():
    def __init__(self) -> None:
        self.accumulate = np.array([0 for _ in range(128)])
        self.colnumbers = [set([]) for _ in range(128)]
        self.PEarray = [[PE(y,x,self) for x in range(128)] for y in range(128)]
        self.cycleCount = 0
    
    def fill_sys(self,compcols, matrix, partitionNumber):
        for compcolNumber,compcol in enumerate(compcols):
            for column in compcol:
                rows = column[1]
                for row in rows:
                    truerow = row + 128*partitionNumber
                    s = list(matrix.indices[matrix.indptr[truerow]:matrix.indptr[truerow+1]])
                    colind = matrix.indptr[truerow] + s.index(column[0])
                    self.colnumbers[compcolNumber].add(column[0])
                    d = matrix.data[colind]
                    self.PEarray[row][compcolNumber].loadval(d,column[0])

    def compute_spmv(self,vector):
        count = 0
        continueFlag = True
        while continueFlag:
            continueFlag = False
            for x in range(127,-1,-1):
                # add in new values from the top:
                if count <= x:
                    if bool(self.colnumbers[x]):
                        col = self.colnumbers[x].pop()
                        self.PEarray[0][x].inputTop(vector[col],col,False)
                    else:
                        self.PEarray[0][x].inputTop(0,0,True)
                for y in range(127,-1,-1):
                    continueFlag |= not self.PEarray[y][x].cycle()

        return self.accumulate
    
    def drain_sys(self):
        self.accumulate = np.array([0 for _ in range(128)])
        self.colnumbers = [set([]) for _ in range(128)]
        arr = self.PEarray
        for x in arr:
            for y in x:
                y.reset()
    
    def print_sys(self):
        arr = self.PEarray
        print("-----beginarray------")
        for x in arr:
            for y in x:
                print(y, end="")
            print()
        print("-----endarray------")
    
    def print_sys_col(self):
        arr = self.PEarray
        print("-----beginarray------")
        for x in arr:
            for y in x:
                print(y.column, end="")
            print()
        print("-----endarray------")



def compressAndCompute(csrm, vector):
    partitionedrows = []
    length = len(csrm.indptr)
    for x in range(math.ceil(length/129)):
        partitionedrows.append(csrm.indptr[x*128:min(length,(x+1)*128+1)])

    
    systolicarray = SystolicArray()
    
    final = []
    

    for (partitionNumber,partition) in enumerate(partitionedrows):
        compressedcols = []

        # Each partition is a set of row pointers that is 129 values long, so encompasses 128 rows.
        # The goal is to go through all columns and potentially compress them.
        # compressedcols: A list of however many partitions there are. Each value in the list is a list 
        # containing compressed columns for a given partition
        # leftovers: values that were rejected from the initial merge
        
        # FOR EACH PARTITION: 
        # create a dictionary of column values, where each holds a set of rows that they contain. 
        # Then, look at the first column in the dictionary, and check all the other columns to see if they fit into it. 
        # If they do not share rows within a given tolerance, then combine a set of columns and delete them from the dictionary
        # Add this compressed column, as well as the leftover row, column values to their respective lists.
        # Repeat this process until the entire dictionary is empty, then repeat for leftover rows.
        # repeat from the start for each partition.
        
        
        #NOTE: datapoints can be extracted from partitions later. Focus on rows for now.
        coldict = {}
        for row in range(0, len(partition)-1):
            colnumbers = csrm.indices[partition[row]: partition[row+1]]
            for x in colnumbers:
                if x in coldict:
                    coldict[x].add((row))
                else:
                    coldict[x] = set([(row)])
        
        while bool(coldict):
            keylist = list(coldict.keys()) # list of all column values in this partition
            baseSet = coldict[keylist[0]].copy() # the rows currently occupied in base.
            
            whichRow = [(keylist[0], coldict[keylist[0]])] # allows us to figure out which rows come from each base.
            whichRowRemainder = []
            coldict.pop(keylist[0])
            keylist = keylist[1:]
            
            intersectedRows = set() # the max intersections for a specific row is implicitly 1.
            totalintersections = 0
            maxintersections = 14
            for x in keylist:
                intersection = baseSet & coldict[x]
                # if there is a row being intersected twice, or the addition would cause us to go over the max intersections for a grouped row:
                if not baseSet.symmetric_difference(coldict[x]) or intersection & intersectedRows or totalintersections + len(intersection) > maxintersections:
                    continue
                    
                
                totalintersections += len(intersection) # increment the total number of intersections
                addedValues = coldict[x] - intersection# row values added this iteration
                intersectedRows |= intersection # add the intersected rows to the variable.
                
                baseSet |= addedValues # add values to baseset
                
                if addedValues:
                    whichRow.append((x, addedValues)) # append the column value, and set of row values added in this iteration
                if intersection:
                    compressedcols.append([(x, intersection)])
                        
                
                coldict.pop(x) # remove the value of x after it has been compressed and split into whichRow and remainders
            compressedcols.append(whichRow)
            compressedcols.append(whichRowRemainder)
        
        
        intermediate = {}
        for x in range(len(compressedcols)):
            for y in compressedcols[x]:
                if x in intermediate:
                    intermediate[x] |= y[1]
                else:
                    intermediate[x] = y[1].copy()
        
        
        #second round of compression
        while bool(intermediate):
            keylist = list(intermediate.keys()) # list of all column values in this partition
            originalkey = keylist[0]
            baseSet = intermediate[keylist[0]].copy() # the rows currently occupied in base.
            
            intermediate.pop(keylist[0])
            keylist = keylist[1:]
                        
            for x in keylist:
                intersection = baseSet & intermediate[x]
                
                # if there is any overlap at all, just skip it.
                if intersection:
                    continue
                
                addedValues = intermediate[x] - intersection# row values added this iteration

                                
                baseSet |= addedValues # add values to baseset
                
                
                
                compressedcols[originalkey] += compressedcols[x]
                
                compressedcols[x] = []

                                                
                intermediate.pop(x)
        # print(compressedcols) # array of arrays filled with (col, row set). This compressed columns in each array, and needs to be unpacked
        # print(intermediate) # second round of compression needs to be done to these leftovers, and then they can be added as additional compressed columns
        
        compressedcols = list(filter(lambda i: i != [], compressedcols))
        
        
        #output is a list of (partition row #, compressed matrix) pairs. The partition number allows us to 
        # figure out the true value of the rows in case we need it, and each column number is preserved anyways.
        # This output will then be partitioned further into compressed columns of 128 each, then fed into the systolic array.
        partial = np.array([0 for _ in range(128)])
        for compressedcolumns in range(0,len(compressedcols),128):
            systolicarray.fill_sys(compressedcols[compressedcolumns*128: min(len(compressedcols), (compressedcolumns + 1) * 128)],csrm,partitionNumber)
            partial += systolicarray.compute_spmv(vector=vector)
        
        print(partial)
        
        final += partial[0:len(vector)- partitionNumber*128].tolist()
    
        systolicarray.drain_sys()
    
    return final
                
                
                
    
    


# Example usage
rows, cols = 200,200
gen = np.random.default_rng()
mask = gen.random((rows, cols)) < 0.1  # Example sparse matrix
vec = gen.integers(low=1, high = 10, size = cols)
data = gen.integers(low=1,high=10,size=(rows, cols)) * mask  # Example sparse matrix
csr_sparse_matrix = csr_matrix(data)
systolic_array_length = 100  # Example length of systolic array
collision_threshold = 0.2  # Example collision threshold

print(csr_sparse_matrix.toarray())
x = compressAndCompute(csr_sparse_matrix,vec)
print(x)

true = np.matmul(csr_sparse_matrix.toarray(),vec.T)
print(true)

print(np.equal(x,true))
