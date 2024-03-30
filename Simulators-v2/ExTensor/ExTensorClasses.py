from collections import deque
import time
class CSFNode:
    def __init__(self, value):
        self.value = value
        self.children = {}
        self.used = False
    
    def __str__(self) -> str:
        print_csf_tree(self)
        return ""

def coo_to_csf(coo_matrix, llbxsize, llbysize, pexsize, peysize, swap):
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
        # NOTICE: b comes before a because THE LLB TILES ARE COLUMN MAJOR, NOT ROW MAJOR. BECAUSE THE MULTIPLICATION IS B-STATIONARY, THE MATRIX MUST BE COL MAJOR FOR THE FIRST TWO COORDINATES
        # WHEN DOING A-STATIONARY, WE WANT BOTH AS ROW-MAJOR, AND WHEN OUTPUT STATIONARY WE WANT A AS ROW MAJOR, AND B AS COL MAJOR
        for coord in [b, a, c, d, e, f]:
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
        
        

class DRAMIntersector:
    def __init__(self, skipto) -> None:
        self.streamOne = []
        self.streamTwo = []
        self.nodeOne = None
        self.nodeTwo = None
        self.LLBIntersector = None
        self.skipto = skipto
        self.endFlag = False
        self.nodeTwoIndices = []
    
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()

        
    
    def setNext(self, llbi: 'LLBIntersector') -> None:
        # The LLBIntersector takes in values of (CSF1, CSF2), where we intersect the children of the CSF2 node with every child of the CSF1 node
        # Due to this, there is a lengthy buffer that the LLBIntersector draws its coordinate values from.
        self.LLBIntersector = llbi
    
    def input(self, s1: CSFNode, s2: CSFNode):
        self.nodeOne = s1
        self.nodeTwo = s2
        self.nodeCounter = 0
        self.nodeTwoIndices = sorted(list(s2.children.keys()))

    def cycle(self) -> None:
        if self.endFlag:
            return
        
        # If one of the streams runs out, we reset node 1 (A) and iterate to the next child of node 2 (B) (output B stationary)
        if not (self.streamOne and self.streamTwo):
            if self.nodeCounter >= len(self.nodeTwo.children):
                self.LLBIntersector.load(None, None, -1, -1)
                self.endFlag = True
                return
            self.streamOne = sorted(list(self.nodeOne.children.keys()),reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter]].children.keys()), reverse=True)
            self.nodeCounter += 1
        
        if self.skipto:
            # Iterate until you find either the end of one of the streams or they equal each other (this should be one popping until it catches up)
            while self.streamOne and self.streamTwo and self.streamOne[-1] != self.streamTwo[-1]:
                if self.streamOne[-1] > self.streamTwo[-1]:
                    self.streamTwo.pop()
                else:
                    self.streamOne.pop()
            
            # if EOF, just go to the next cycle to load the next values.
            if not (self.streamOne and self.streamTwo):
                return
            
            # load node One's children one level down, and node Two's child two levels down. This is because both sides are in row major form (HAVE TO CHECK CSF DEFINITION TO MAKE SURE OF THIS)
            # In LLBIntersector, we take in nodeOne's Child (a "Node"), and intersect each of its children with the children of node Two's child two levels down (a "Scalar")
            node2Col = self.nodeTwoIndices[self.nodeCounter-1]
            node2 = self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter-1]].children[self.streamTwo.pop()]
            node1List = self.nodeOne.children[self.streamOne.pop()].children.items()
            for node in node1List:
                #print("newSet")
                #print("i1")
                #print_csf_tree(node[1])
                #print("i2")
                #print_csf_tree(node2)
                self.LLBIntersector.load(node[1], node2, node[0], node2Col )

        else:
            #Similar logic to the skipTo side above, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                node2Col = self.nodeTwoIndices[self.nodeCounter-1]
                node2 = self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter-1]].children[self.streamTwo.pop()]
                node1List = self.nodeOne.children[self.streamOne.pop()].children.items()
                for node in node1List:
                    self.LLBIntersector.load(node[1], node2, node[0], node2Col)            
            elif self.streamOne[-1] > self.streamTwo[-1]:
                self.streamTwo.pop()
            else:
                self.streamOne.pop()
    

class LLBIntersector:
    def __init__(self, skipto) -> None:
        self.streamOne = []
        self.streamTwo = []
        self.nodeOne = None
        self.nodeTwo = None
        self.PEArray = None
        self.skipto = skipto
        self.endFlag = False
        self.nodeOneIndices = []
        self.inputBuffer = deque([])
        self.inputCoordBuffer = deque([])
        self.inputFlag = False
        self.nodeCounter = 0
        self.cycleFlag = False

    def setNext(self,PEArray) -> None:
        self.PEArray = PEArray
        
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def load(self, s1: CSFNode|None, s2: CSFNode|None, LLBrow, LLBcol) -> None:
        self.inputBuffer.append((s1,s2))
        self.inputCoordBuffer.append((LLBrow,LLBcol))
        self.inputFlag = True
    
    def cycle(self) -> None:
        if self.endFlag or not self.inputFlag:
            return
        
        # If one of the streams runs out, we reset node 1 (A) and iterate to the next child of node 2 (B) (output B stationary)
        if not (self.streamOne and self.streamTwo):
            if not (self.nodeOne or self.nodeTwo) or self.nodeCounter >= len(self.nodeOne.children):
                if not self.inputBuffer:
                    self.inputFlag = False
                    return
                newNodes = self.inputBuffer.popleft()
                self.row, self.col  = self.inputCoordBuffer.popleft()
                if newNodes == (None,None):
                    self.endFlag = True
                    self.PEArray.load(None,None,None,None,None,None)
                    return
                self.nodeOne = newNodes[0]
                self.nodeTwo = newNodes[1]
                self.nodeCounter = 0
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                return
            self.streamOne = sorted(list(self.nodeOne.children[self.nodeOneIndices[self.nodeCounter]].children.keys()), reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children.keys()),reverse=True)
            self.nodeCounter += 1
        if self.skipto:
            # Iterate until you find either the end of one of the streams or they equal each other (this should be one popping until it catches up)
            while self.streamOne and self.streamTwo and self.streamOne[-1] != self.streamTwo[-1]:
                if self.streamOne[-1] > self.streamTwo[-1]:
                    self.streamTwo.pop()
                else:
                    self.streamOne.pop()

            # if EOF, just go to the next cycle to load the next values.
            if not (self.streamOne and self.streamTwo):
                return
            

            node1Row = self.nodeOneIndices[self.nodeCounter-1]
            node = self.nodeOne.children[self.nodeOneIndices[self.nodeCounter-1]].children[self.streamOne.pop()]
            node2List = self.nodeTwo.children[self.streamTwo.pop()].children.items()
            for node2 in node2List:
                #print("newSet")
                #print("i1")
                #print_csf_tree(node)
                #print("i2")
                #print_csf_tree(node2[1])
                self.PEArray.load(node, node2[1], self.row, self.col, node1Row, node2[0])

        else:
            #Similar logic to the skipTo side above, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                node1Row = self.nodeOneIndices[self.nodeCounter-1]
                node = self.nodeOne.children[self.nodeOneIndices[self.nodeCounter-1]].children[self.streamOne.pop()]
                node2List = self.nodeTwo.children[self.streamTwo.pop()].children.Items()
                for node2 in node2List:
                    self.PEArray.load(node, node2[1], self.row, self.col, node1Row, node2[0])
            elif self.streamOne[-1] > self.streamTwo[-1]:
                self.streamTwo.pop()
            else:
                self.streamOne.pop()
    

class PEArray:
    def __init__(self, npe:int) -> None:
        self.i = None
        self.nexti = None
        self.endFlag = True
        self.inputBuffer = deque([])
        self.currentPE = 0
        self.numPEs = npe
        self.endFlag = False
        self.inputFlag = False
        self.PEs = []
    
    def setNext(self, PEs: list['PEIntersector']):
        self.PEs = PEs
        
        
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def load(self, s1: CSFNode|None, s2: CSFNode|None, LLBrow, LLBCol, PERow, PECol) -> None:
        self.inputBuffer.append((s1,s2, LLBrow, LLBCol, PERow, PECol))
        self.inputFlag = True
    
    def cycle(self) -> None:
        if self.endFlag or not self.inputFlag:
            return
        if not self.inputBuffer:
            self.inputFlag = False
            return
        i = self.inputBuffer.popleft()
        if i[0] == None or i[1] == None:
            self.endFlag = True
            for x in self.PEs:
                x.load((None,None,None,None))
        self.PEs[self.currentPE].load(i)
        self.currentPE += 1
        self.currentPE %= self.numPEs

    

class PEIntersector:  
    def __init__(self, skipto, llbTileSize, PETileSize, num) -> None:
        self.streamOne = []
        self.streamTwo = []
        self.nodeOne = None
        self.nodeTwo = None
        self.output = []
        self.skipto = skipto
        self.endFlag = False
        self.inputFlag = False
        self.nodeOneIndices = []
        self.nodeCounter1 = 0
        self.nodeCounter2 = 0
        self.inputBuffer = deque([])
        self.inputCoordBuffer = deque([])
        self.llbTileSize = llbTileSize
        self.peTileSize = PETileSize
        self.id = num
    
    def load(self, input: tuple) -> None:
        # put the offsets for the coordinates here. The final coords can be added to this to get the x and y of the output.
        if input[0] != None:
            #print("PRINTING INPUT TO: ", self.id)
            #print("I1: ")
            #print_csf_tree(input[0])
            #print("I2:", )
            #print_csf_tree(input[1])
            self.inputBuffer.append((input[0],input[1]))
            self.inputCoordBuffer.append((self.llbTileSize * input[2] + self.peTileSize * input[4], self.llbTileSize * input[3] + self.peTileSize * input[5]))
            self.inputFlag = True
        else:
            self.inputBuffer.append((None,None))
            self.inputCoordBuffer.append((0,0))
            self.inputFlag = True
            
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def cycle(self) -> None:   
        if self.endFlag or not self.inputFlag:
            return
        
        # If one of the streams runs out, we reset node 1 (A) and iterate to the next child of node 2 (B) (output B stationary)
        if not (self.streamOne and self.streamTwo):
            if self.nodeOne and self.nodeTwo and self.nodeCounter2 >= len(self.nodeTwo.children) and self.nodeCounter1 < len(self.nodeOne.children):
                self.nodeCounter1 += 1
                self.nodeCounter2 = 0
            elif not (self.nodeOne or self.nodeTwo) or self.nodeCounter2 >= len(self.nodeTwo.children):
                if not self.inputBuffer:
                    self.inputFlag = False
                    return
                newNodes = self.inputBuffer.popleft()
                self.row, self.col  = self.inputCoordBuffer.popleft()
                if newNodes == (None,None):
                    self.endFlag = True
                    return
                self.nodeOne = newNodes[0]
                self.nodeTwo = newNodes[1]
                self.nodeCounter2 = 0
                self.nodeCounter1 = 1
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                self.nodeTwoIndices = sorted(list(self.nodeTwo.children.keys()))
                return
            self.nodeCounter2 += 1
            self.streamOne = sorted(list(self.nodeOne.children[self.nodeOneIndices[self.nodeCounter1-1]].children.keys()),reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter2-1]].children.keys()),reverse=True)
            return
        
        if self.skipto:
            # Iterate until you find either the end of one of the streams or they equal each other (this should be one popping until it catches up)
            while self.streamOne and self.streamTwo and self.streamOne[-1] != self.streamTwo[-1]:
                if self.streamOne[-1] > self.streamTwo[-1]:
                    self.streamTwo.pop()
                else:
                    self.streamOne.pop()

            # if EOF, just go to the next cycle to load the next values.
            if not (self.streamOne and self.streamTwo):
                return
            

            node1Row = self.nodeOneIndices[self.nodeCounter1-1]
            node2Col = self.nodeTwoIndices[self.nodeCounter2-1]
            node1 = self.nodeOne.children[self.nodeOneIndices[self.nodeCounter1-1]].children[self.streamOne.pop()]
            node2 = self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter2-1]].children[self.streamTwo.pop()]
            self.output.append((node1.value * node2.value, self.row + node1Row, self.col + node2Col))

        else:
            #Similar logic to the skipTo side above, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                node1Row = self.nodeOneIndices[self.nodeCounter-1]
                node = self.nodeOne.children[self.nodeOneIndices[self.nodeCounter-1]].children[self.streamOne.pop()]
                node2List = self.nodeTwo.children[self.streamTwo.pop()].children.Items()
                for node2 in node2List:
                    self.PEArray.load(node, node2[1], self.row, self.col, node1Row, node2[0])
            elif self.streamOne[-1] > self.streamTwo[-1]:
                self.streamTwo.pop()
            else:
                self.streamOne.pop()