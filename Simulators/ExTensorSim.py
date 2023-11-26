from typing import Callable
import numpy as np
import scipy.io as sio
from scipy.sparse import csc_array
from scipy.sparse import csr_array
from scipy.sparse import coo_matrix


sparseMat = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')

class CSFNode:
    def __init__(self, value):
        self.value = value
        self.children = {}
    
    def __str__(self) -> str:
        print_csf_tree(self)
        return "Root"

def coo_to_csf(coo_matrix, LLBxsize, LLBysize, PExsize, PEysize,swap):
    root = CSFNode(None)

    for i, j, data in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
        current_node = root
        
        a = i // LLBxsize
        b = j // LLBysize
        
        c = (i - (a * LLBxsize)) // PExsize
        d = (j - (b * LLBysize)) // PExsize
        
        
        e = (i - (a * LLBxsize) - (c * PExsize))
        f = (j - (b * LLBysize) - (d * PEysize))
        
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
"""

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
    
        

    

class Intersect:  
    def reset(self) -> None:
        self.A = None
        self.B = None
        self.skipTo = False
          
    def set(self,scanner1: Scanner, scanner2: Scanner, skipTo: bool ) -> None:
        self.A = scanner1
        self.B = scanner2
        self.skipTo = skipTo
        
    
    def getScanner1(self) -> Scanner:
        return self.A
    
    def getScanner2(self) -> Scanner:
        return self.B

    def Intersect(self):
        #NOTE: SkipTo architecture is not directly implemented, and only matters for cycle counting.
        output = []
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
                output.append(self.A.iterate())
        return output

                

# Node 1 Stationary.
def nodeStationary(node1: CSFNode, node2: CSFNode, level, reverse):
    intersector = Intersect()
    scanner1 = Scanner()
    scanner2 = Scanner()
    total = []
    for i in node1.children:
        value = node1.children[i]
        scanner1.loadData(list(value.children.keys()),level+2)
        
        scanner2.loadData(list(node2.children.keys()),level+1)
        
        intersector.set(scanner1,scanner2, True)
        
        for k in intersector.Intersect():
            for j in node2.children[k].children:
                if reverse:
                    total.append((node2.children[k].children[j],value.children[k],i,j))
                else:
                    total.append((value.children[k],node2.children[k].children[j],i,j))
            
    return total

def outputStationary(node1: CSFNode, node2: CSFNode, level):
    intersector = Intersect()
    scanner1 = Scanner()
    scanner2 = Scanner()
    total = []
    for i in node1.children:
        for j in node2.children:
            value1 = node1.children[i]
            value2 = node2.children[j]
            scanner1.loadData(list(value1.children.keys()),level + 2)
            scanner2.loadData(list(value2.children.keys()),level + 2)
            intersector.set(scanner1,scanner2, True)
            
            for k in intersector.Intersect():
                total.append((value1.children[k],value2.children[k],i,j))
    return total

"""
gen = np.random.default_rng()
data1 = gen.integers(0,10,5)
row1 = gen.integers(0,5,5)
col1 = gen.integers(0,5,5)

data2 = gen.integers(0,10,5)
row2 = gen.integers(0,5,5)
col2 = gen.integers(0,5,5)
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
output1 = nodeStationary(input1,input2,0,False)
for x in output1:
    output2 = nodeStationary(x[0],x[1],2,True)
    for y in output2:
        output3 = outputStationary(y[0],y[1],4)
        for z in output3:
            datalist.append(z[0].value * z[1].value)
            rowlist.append(x[2]*3 + y[2]*2 + z[3])
            collist.append(x[3]*3 + y[3]*2 + z[2])

print(rowlist)
print(collist)
print(len(datalist))
out = coo_matrix((datalist,(rowlist,collist)), shape=(5,5))
print("out\n",out.toarray())
print("matmul\n",np.matmul(i1.toarray(),i2.toarray()))
"""

datalist = []
rowlist = []
collist = []

output1 = nodeStationary(coo_to_csf(sparseMat, 25, 25, 5, 5,False), coo_to_csf(sparseMat, 25, 25, 5, 5,True),0,False)
for x in output1:
    output2 = nodeStationary(x[0],x[1],2,True)
    for y in output2:
        output3 = outputStationary(y[0],y[1],4)
        for z in output3:
            datalist.append(z[0].value * z[1].value)
            # switch up z values since they were swapped in the original array.
            rowlist.append(x[2]*25 + y[2]*5 + z[3])
            collist.append(x[3]*25 + y[3]*5 + z[2])
        
    #datalist.append(x[0].value * x[1].value)
    #rowlist.append(x[2])
    #collist.append(x[3])
out = coo_matrix((datalist,(rowlist,collist)), shape=(496,496))
print(out.toarray())
print(np.matmul(sparseMat.toarray(),sparseMat.toarray()))
print(np.array_equal(out.toarray() == 0, np.matmul(sparseMat.toarray(),sparseMat.toarray()) == 0 ))
print(np.allclose(out.toarray(),np.matmul(sparseMat.toarray(),sparseMat.toarray()),rtol=0.000001))
