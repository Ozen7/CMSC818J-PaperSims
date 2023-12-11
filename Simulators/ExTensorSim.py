import numpy as np
import scipy.io as sio
from scipy.sparse import coo_matrix
import bisect


sparseMat = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')
global numCycles

class CSFNode:
    def __init__(self, value):
        self.value = value
        self.children = {}
        self.used = False
    
    def __str__(self) -> str:
        print_csf_tree(self)
        return "Root"

def coo_to_csf(coo_matrix, llbxsize, llbysize, pexsize, peysize,swap):
    root = CSFNode(None)

    for i, j, data in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
        current_node = root
        
        a = i // llbxsize
        b = j // llbysize
        
        c = (i - (a * llbxsize)) // pexsize
        d = (j - (b * llbysize)) // pexsize
        
        
        e = (i - (a * llbxsize) - (c * pexsize))
        f = (j - (b * llbysize) - (d * peysize))
        
        if swap:
            temp = e
            e = f
            f = temp

        # Traverse the tree based on the coordinates
        for coord in [a, b, c, d, e, f]:
            if coord not in current_node.children:
                current_node.children[coord] = CSFNode(None)
            current_node = current_node.children[coord]

        # Set the leaf node value to the data value
        if current_node.value != None:
            current_node.value += data
        else:
            current_node.value = data

    return root

# Print the CSF representation (basic representation)
def print_csf_tree(node, depth=0):
    if node.value is not None:
        print(f"{'  ' * depth}Leaf: {node.value}")
    for coord, child in node.children.items():
        print(f"{'  ' * depth}Coordinate {coord}:")
        print_csf_tree(child, depth + 1)
"""   
# Create a sample COO matrix
data = [1, 2, 3, 4, 5,10]
row = [0, 1, 1, 2, 1,5]
col = [0, 0, 1, 2, 2,5]

coo_matrix_example = coo_matrix((data, (row, col)), shape=(6, 6))

# Convert COO to CSF
csf_representation = coo_to_csf(coo_matrix_example)

#print_csf_tree(csf_representation)


class Scanner:
    
    metadata = []
    level = 0
    
    # Metadata is simply an array of coordinates. By repeatedly intersecting two coordinate streams, we can obtain the 
    # final pairing of operations necessary to get the output.
    def loadData(self, metadata, level) -> None:
        self.metadata = sorted(metadata) 
        self.metadata.reverse() # to make pop() faster
        self.level = level    
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
"""
        

    

class Intersect:  
    def __init__(self) -> None:
        self.A = []
        self.B = []
        self.skipto = False
          
    def set(self,scanner1: list, scanner2: list, skipto: bool, peNumber: int) -> None:
        global numCycles
        numCycles[peNumber] += 1
        self.peNumber = peNumber
        self.A = sorted(scanner1)
        self.B = sorted(scanner2)
        self.skipto = skipto
        
    def skipToA(self,target: int):
        index = bisect.bisect_left(self.A, target)
        if index < len(self.A):
            self.A = self.A[index:]
        else:
            self.A = []
    
    def skipToB(self,target: int):
        index = bisect.bisect_left(self.B, target)
        if index < len(self.B):
            self.B = self.B[index:]
        else:
            self.B = []
        
    
    def getScanner1(self) -> list:
        return self.A
    
    def getScanner2(self) -> list:
        return self.B

    def Intersect(self):
        global numCycles
        output = []
        if self.skipto:
            while True:
                numCycles[self.peNumber] += 1
                if self.A == [] or self.B == []:
                    self.A = []
                    self.B = []
                    break
                elif self.A[0] == self.B[0]:
                    self.B.pop(0)
                    output.append(self.A.pop(0))
                else:
                    temp = self.A[0]
                    self.skipToA(self.B[0])
                    self.skipToB(temp)
        else:
            while True:
                numCycles[self.peNumber] += 1
                if self.A == [] or self.B == []:
                    self.A = []
                    self.B = []
                    break
                elif self.A[0] < self.B[0]:
                    self.A.pop(0);
                elif self.A[0] > self.B[0]:
                    self.B.pop(0);
                else:
                    self.B.pop(0)
                    output.append(self.A.pop(0))
        return output

                

# Node 1 Stationary.
def nodeStationary(node1: CSFNode, node2: CSFNode, reverse, peNumber):
    node1.used = True
    node2.used = True
    intersector = Intersect()
    total = []
    for i in node1.children:
        value = node1.children[i]
        value.used = True

        intersector.set(list(value.children.keys()),list(node2.children.keys()), skipto=True, peNumber=peNumber)
        
        for k in intersector.Intersect():
            for j in node2.children[k].children:
                node2.children[k].used = True
                # add in reversing so I don't have to deal with the headache of indices swapping 
                # when I want to mix matrix A and B stationary
                if reverse:
                    total.append((node2.children[k].children[j],value.children[k],i,j))
                else:
                    total.append((value.children[k],node2.children[k].children[j],i,j))
            
    return total

def outputStationary(node1: CSFNode, node2: CSFNode,peNumber):
    node1.used = True
    node2.used = True
    intersector = Intersect()
    total = []
    for i in node1.children:
        for j in node2.children:
            value1 = node1.children[i]
            value2 = node2.children[j]
            value1.used = True
            value2.used = True
            intersector.set(list(value1.children.keys()),list(value2.children.keys()), skipto=True, peNumber=peNumber)
            for k in intersector.Intersect():
                total.append((value1.children[k],value2.children[k],i,j))
    return total

"""
gen = np.random.default_rng()
data1 = gen.integers(0,10,10)
row1 = gen.integers(0,5,10)
col1 = gen.integers(0,5,10)

data2 = gen.integers(0,10,10)
row2 = gen.integers(0,5,10)
col2 = gen.integers(0,5,10)
i1 = coo_matrix((data1, (row1, col1)), shape=(5, 5))
i2 = coo_matrix((data2, (row2, col2)), shape=(5, 5))

input1 = coo_to_csf(i1,3,3,2,2,False)
input2 = coo_to_csf(i2,3,3,2,2,True)

#print("i1")
#print_csf_tree(input1)
#print("i2")
#print_csf_tree(input2)

print(i2.toarray())
print(i1.toarray())


datalist = []
rowlist = []
collist = []
output1 = nodeStationary(input1,input2,False)
for x in output1:
    output2 = nodeStationary(x[0],x[1],True)
    for y in output2:
        output3 = outputStationary(y[0],y[1])
        for z in output3:
            numCycles += 1 # assuming multiplication takes one cycle
            datalist.append(z[0].value * z[1].value)
            rowlist.append(x[2]*3 + y[2]*2 + z[3])
            collist.append(x[3]*3 + y[3]*2 + z[2])

out = coo_matrix((datalist,(rowlist,collist)), shape=(5,5))
print("out\n",out.toarray())
print("matmul\n",np.matmul(i1.toarray(),i2.toarray()))

"""

def count_csf_tree(node):
    intsLoaded = 0
    if node.used:
        node.used = False
        intsLoaded += 1
        for child in node.children.values():
            intsLoaded += count_csf_tree(child)
    return intsLoaded

def run_extensor(matrix1, matrix2, l1tile, l2tile):
    datalist = []
    rowlist = []
    collist = []
    global numCycles
    numCycles = [0 for _ in range(12)] # 1 DRAM, 1 LLB, 10 PEs
    lastIterCycles = [0 for _ in range(10)]
    datacount = 0
    totalcycles = 0
    coo_m1 = coo_matrix(matrix1)
    coo_m2 = coo_matrix(matrix2)
    input1 = coo_to_csf(coo_m1, l1tile, l1tile, l2tile, l2tile,False)
    input2 = coo_to_csf(coo_m2, l1tile, l1tile, l2tile, l2tile,True)
    output1 = nodeStationary(input1, input2,False,0)
    for x in output1:
        datacount += 2
        output2 = nodeStationary(x[0],x[1],True,1)
        count = 0
        for y in output2:
            lastIterCycles = [0 for _ in range(10)]
            output3 = outputStationary(y[0],y[1],2 + np.argmin(numCycles[2:]))
            count += 1
            for z in output3:
                lastIterCycles[ count%10 ] += 1
                datalist.append(z[0].value * z[1].value)
                # switch up z values since they were swapped in the original array.
                rowlist.append(x[2]*l1tile + y[2]*l2tile + z[3])
                collist.append(x[3]*l1tile + y[3]*l2tile + z[2])
        
        datacount += count_csf_tree(x[0])
        datacount += count_csf_tree(x[1])
        
        if max(numCycles[2:12]) > numCycles[1]:
            totalcycles += max(numCycles[2:12])
        else:
            totalcycles += numCycles[1] + max(lastIterCycles)
            
        numCycles = [0 for _ in range(12)]
        lastIterCycles = [0 for _ in range(10)]
            
        
        #datalist.append(x[0].value * x[1].value)
        #rowlist.append(x[2])
        #collist.append(x[3])
    out = coo_matrix((datalist,(rowlist,collist)), shape=(matrix1.shape[0],matrix2.shape[1]))
    # Calculate the number of integers loaded from memory by going through the CSFNode tree and counting the number of used nodes.

    print("Cycles:", numCycles)
    print("Total Cycles (parallel)", totalcycles)
    print("Integers Loaded From Memory:", datacount)
    if matrix1.size >= 10000 or matrix2.size >= 10000:
        print("matrices too large to verify")
    else:
        trueval = matrix1 @ matrix2
        print("Verify that sparse multiplication is correct: ", np.allclose(out.toarray(),trueval.toarray(),rtol=0.000001))
 
 
'''
gen = np.random.default_rng()
data1 = gen.integers(1,10,100)
row1 = gen.integers(0,1000,100)
col1 = gen.integers(0,1000,100)

data2 = gen.integers(1,10,100)
row2 = gen.integers(0,1000,100)
col2 = gen.integers(0,1000,100)
i1 = coo_matrix((data1, (row1, col1)), shape=(1000, 1000))
i2 = coo_matrix((data2, (row2, col2)), shape=(1000, 1000))
print(i1.toarray())
print(i2.toarray())
run_extensor(i1,i2,250,125)
'''