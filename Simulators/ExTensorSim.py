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
        return coo_to_csf(self)

def coo_to_csf(coo_matrix):
    root = CSFNode(None)

    for i, j, data in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
        current_node = root

        # Traverse the tree based on the coordinates
        for coord in (i, j):
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
    parent_node = 0 
    
    # Metadata is simply an array of coordinates. By repeatedly intersecting two coordinate streams, we can obtain the 
    # final pairing of operations necessary to get the output.
    def loadData(self, metadata, level, parent_node) -> None:
        self.metadata = sorted(metadata) 
        self.metadata.reverse() # to make pop() faster
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

                

      
def outputStationary(node1: CSFNode,node2: CSFNode,level,node2p):
    intersector = Intersect()
    scanner1 = Scanner()
    scanner2 = Scanner()
    total = []
    for i in node1.children:
        value = node1.children[i]
        scanner1.loadData(list(value.children.keys()),level+2,i)
        
        scanner2.loadData(list(node2.children.keys()),level+1,node2p)
        
        intersector.set(scanner1,scanner2, True)
        
        for k in intersector.Intersect():
            for j in node2.children[k].children:
                total.append((value.children[k],node2.children[k].children[j],i,j))
            
    return total

"""
gen = np.random.default_rng()
data1 = gen.integers(0,100,5)
row1 = gen.integers(0,3,5)
col1 = gen.integers(0,3,5)

data2 = gen.integers(0,100,5)
row2 = gen.integers(0,3,5)
col2 = gen.integers(0,3,5)
i1 = coo_matrix((data1, (row1, col1)), shape=(3, 3))
i2 = coo_matrix((data2, (row2, col2)), shape=(3, 3))

input1 = coo_to_csf(i1)
input2 = coo_to_csf(i2)

#print("i1")
#print_csf_tree(input1)
#print("i2")
#print_csf_tree(input2)

print(i2.toarray())
print(i1.toarray())


datalist = []
rowlist = []
collist = []
#inter = outputStationary(input1,input2,0,0)
#for x in inter:
#    print(x[0].value, x[1].value,x[2],x[3])

    
#print(len(datalist))
#out = coo_matrix((datalist,(rowlist,collist)), shape=(3,3))
#print("out\n",out.toarray())
#print("matmul\n",np.matmul(i1.toarray(),i2.toarray()))
"""

datalist = []
rowlist = []
collist = []
for x in outputStationary(coo_to_csf(sparseMat),coo_to_csf(sparseMat),0,0):
    datalist.append(x[0].value * x[1].value)
    rowlist.append(x[2])
    collist.append(x[3])

out = coo_matrix((datalist,(rowlist,collist)), shape=(496,496))
print(out.toarray())
print(np.matmul(sparseMat.toarray(),sparseMat.toarray()))

print(np.allclose(out.toarray(),np.matmul(sparseMat.toarray(),sparseMat.toarray()),rtol=0.000001))
 