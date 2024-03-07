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
        
        

class DRAMIntersector: # NEEDS TESTING
    def __init__(self, skipto) -> None:
        self.streamOne = []
        self.streamTwo = []
        self.nodeOne = None
        self.nodeTwo = None
        self.LLBIntersector = None
        self.skipto = skipto
    
    def setNext(self, llbi: 'LLBIntersector') -> None:
        # The LLBIntersector takes in values of (CSF1, CSF2), where we intersect the children of the CSF2 node with every child of the CSF1 node
        # Due to this, there is a lengthy buffer that the LLBIntersector draws its coordinate values from.
        self.LLBIntersector = llbi
    
    def input(self, s1: CSFNode, s2: CSFNode):
        self.nodeOne = s1
        self.nodeTwo = s2
        self.nodeCounter = 0

    def cycle(self) -> None:
        if self.endFlag:
            return
        
        # If one of the streams runs out, we reset node 1 (A) and iterate to the next child of node 2 (B) (output B stationary)
        if not (self.streamOne and self.streamTwo):
            if self.nodeCounter > len(self.nodeTwo.children):
                self.LLBIntersector.load(None, None)
                self.endFlag = True
                return
            self.streamOne = self.nodeOne.children.keys().sort(reverse = True)
            self.streamTwo = self.nodeTwo[self.nodeCounter].children.keys().sort(reverse = True)
            nodeCounter += 1
        
        if self.skipto:
            # Iterate until you find either the end of one of the streams or they equal each other (this should be one popping until it catches up)
            while self.streamOne and self.streamTwo and self.streamOne[-1] != self.streamTwo[-1]:
                if self.streamOne[-1] > self.streamTwo[-1]:
                    self.streamTwo.pop()
                else:
                    self.streamOne.pop()
            
            # if EOF, just go to the next cycle to load the next values.
            if self.streamOne or self.streamTwo:
                return
            
            # load node One's children one level down, and node Two's child two levels down. This is because both sides are in row major form (HAVE TO CHECK CSF DEFINITION TO MAKE SURE OF THIS)
            # In LLBIntersector, we take in nodeOne's Child (a "Node"), and intersect each of its children with the children of node Two's child two levels down (a "Scalar")
            self.LLBIntersector.load(self.nodeOne.children[self.streamOne.pop()], self.nodeTwo.children[self.nodeCounter].children[self.streamTwo.pop()])
        else:
            #Similar logic to the skipTo side above, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                self.LLBIntersector.load(self.nodeOne.children[self.streamOne.pop()], self.nodeTwo.children[self.nodeCounter].children[self.streamTwo.pop()])
            elif self.streamOne[-1] > self.streamTwo[-1]:
                self.streamTwo.pop()
            else:
                self.streamOne.pop()
    

class LLBIntersector:
    def __init__(self) -> None:
        self.nodeOne = None
        self.nodeTwo = None

    def setNext(self,PEArray) -> None:
        self.PEArray = PEArray
    
    def inputBuffer(self, i):
        self.nexti = i
    
    def load(self, s1: CSFNode|None, s2: CSFNode|None) -> None:
        self.nodeOne = s1
        self.nodeTwo = s2
    
    def cycle(self) -> None:
        if self.i == None:
            return False
        print("LLBI")
        self.PEArray.inputBuffer(self.i)
        return self.i
    

class PEArray:
    def __init__(self) -> None:
        self.i = None
        self.nexti = None
    
    def inputBuffer(self, i):
        self.nexti = i
    
    def setNext(self, PEs: list):
        self.PEs = PEs
    
    def load(self) -> None:
        self.i = self.nexti
        self.nexti = None
    
    def cycle(self) -> None:
        if self.i == None:
            return False
        for x in self.PEs:
            x.inputBuffer(self.i)
        print("PEA")
        return self.i

    

class PEIntersector:  
    def __init__(self) -> None:
        self.i = None
        self.nexti = None

    def setNext(self,output) -> None:
        self.output = output
    
    def inputBuffer(self, i):
        self.nexti = i
    
    def load(self) -> None:
        self.i = self.nexti
        self.nexti = None
    
    def cycle(self) -> bool:
        if self.i == None:
            return False
        print("PEI")
        return self.i