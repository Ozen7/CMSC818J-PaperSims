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
       
       
class Memory:
    def __init__(self) -> None:
        self.requestQueue = {"LLBI":deque([]), "DRAMI":deque([])}
        self.order = deque([])
    
    def requestMemory(self,bytes,name):
        self.requestQueue[name].append(bytes)
        self.order.append(name)
    
    def cycle(self):
        if not self.order:
            return
        self.requestQueue[self.order[0]][0] -= 1128 # the number of bytes sent per cycle by DRAM assuming 1GhZ clock and 128 gb/s
        if self.requestQueue[self.order[0]][0] < 0:
            self.requestQueue[self.order[0]].popleft()
            self.order.popleft()
        
        

class DRAMIntersector:
    def __init__(self, skipto) -> None:
        self.streamOne = []
        self.streamTwo = []
        self.nodeOne = None
        self.nodeTwo = None
        self.nodeCounter = 0
        self.LLBIntersector = None
        self.memory = None
        self.skipto = skipto
        self.endFlag = False
        self.memoryFlag = False
        self.nodeTwoIndices = []
        self.memoryAccessBytes = 0
        self.memoryWastedCycles = 0
        self.T = 32
        
    def __str__(self) -> str:
        return "DRAM" + " StreamOne: " + str(self.streamOne) + ", StreamTwo: " + str(self.streamTwo) + ", endFlag: " + str(self.endFlag)
    
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()

        
    
    def setNext(self, llbi: 'LLBIntersector', memory: 'Memory') -> None:
        # The LLBIntersector takes in values of (CSF1, CSF2), where we intersect the children of the CSF2 node with every child of the CSF1 node
        # Due to this, there is a lengthy buffer that the LLBIntersector draws its coordinate values from.
        self.LLBIntersector = llbi
        self.memory = memory
    
    def input(self, s1: CSFNode, s2: CSFNode):
        self.nodeOne = s1
        self.nodeTwo = s2
        self.nodeCounter = 0
        self.nodeTwoIndices = sorted(list(s2.children.keys()))
        # load in the LLB metadata!:
        coordCounter = 0
        for x in s1.children:
            coordCounter += 1
            for y in s1.children[x].children:
                coordCounter += 1
                
        for x in s2.children:
            coordCounter += 1
            for y in s1.children[x].children:
                coordCounter += 1
                
        self.MemoryAccess(coordCounter * 2) # two bytes per coordinate!
        # Normally, memory access would be done here! Both nodes are loaded up to the LLB metadata (due to post-intersection fill)
        # Its easier to count when individual rows are 
        
    def MemoryAccess(self, numBytes):
        self.memoryAccessBytes += numBytes
        self.memory.requestMemory(numBytes,"DRAMI")
        self.memoryFlag = True
        
        

    def cycle(self) -> None:
        if self.endFlag:                
            return
        
        # if there are things being pulled from memory, don't cycle.
        if self.memoryFlag:
            self.memoryWastedCycles += 1
            if self.memory.requestQueue["DRAMI"]:
                return
            else:
                self.memoryFlag = False
        
        # Critical note about this part of the accelerator: 
        # In the paper it states that the DRAM intersector simply works as a bootstrapper which does not do any intersection, 
        # and instead sends sets of LLB tiles to the next part (even if they are 0-valued). Normally, you would just match tile (0,1) with 
        # Column 1 in B, (0,2) with Column 2, etc... but since the data is in CSF form, if Column 1 is 0-valued, it simply doesn't exist
        # within the data, meaning you can't just do a 1-1 match. This means that the ONLY way to find proper matches that can be 
        # sent to the LLB tiles is to use intersection logic, since you need to "search" for matching columns for a given tile!
        # If all LLB tiles are nonzero, then this works exactly like a bootstrapper with no intersection logic,
        # but for the sake of avoiding issues with mismatching, there is intersection logic at this level.
        
        # If one of the streams runs out, we reset node 1 (A) and iterate to the next child of node 2 (B) (output B stationary)
        # This is done in a post-intersection fill manner. (Only the arrays themselves are pulled from memory, no other data is necessary)
        # When sending information to the LLB intersector, only the LLB metadata is used as well, meaning the LLB intersector has to pull
        # PE level metadata from memory to do its iteration (also post-intersection fill) 

        if not (self.streamOne and self.streamTwo):
            if self.nodeCounter >= len(self.nodeTwo.children):
                self.LLBIntersector.load(None, None, -1, -1)
                self.endFlag = True
                return
            self.streamOne = sorted(list(self.nodeOne.children.keys()),reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter]].children.keys()), reverse=True)
            self.nodeCounter += 1
        
        if self.skipto:
            counter = 0
            if self.streamOne[-1] < self.streamTwo[-1]:
                while (self.streamOne and self.streamTwo) and self.streamOne[-1] < self.streamTwo[-1] and counter < self.T:
                    self.streamOne.pop()
                    counter += 1
            elif self.streamOne[-1] > self.streamTwo[-1]:
                while (self.streamOne and self.streamTwo) and self.streamOne[-1] > self.streamTwo[-1] and counter < self.T:
                    self.streamTwo.pop()
                    counter += 1
    
            # if EOF or not equal, just go to the next cycle to load the next values.
            if not (self.streamOne and self.streamTwo) or self.streamOne[-1] != self.streamTwo[-1]:
                return
            
            # load node One's children one level down, and node Two's child two levels down. This is because both sides are in row major form (HAVE TO CHECK CSF DEFINITION TO MAKE SURE OF THIS)
            # In LLBIntersector, we take in nodeOne's Child (a "Node"), and intersect each of its children with the children of node Two's child two levels down (a "Scalar")
            node2Col = self.nodeTwoIndices[self.nodeCounter-1]
            node2 = self.nodeTwo.children[self.nodeTwoIndices[self.nodeCounter-1]].children[self.streamTwo.pop()]
            node1List = self.nodeOne.children[self.streamOne.pop()].children.items()
            for node in node1List:
                # The DRAM intersector doesn't send ANYTHING but coordinates to the LLB intersector, which pulls necessary info from DRAM
                # According to the paper, it doesn't have data storage or intersection logic, so the only thing it really can send is coordinate data, unless it loads PE level data 
                # From DRAM when its sending. Either way, I have it setup so that the LLB intersector pulls the data from memory as soon as it loads it in (which runs in the DRAM thread
                # regardless, so no real slow down.)
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
        self.memory = None
        self.skipto = skipto
        self.endFlag = False
        self.nodeOneIndices = []
        self.inputBuffer = deque([])
        self.inputCoordBuffer = deque([])
        self.inputFlag = False
        self.memoryFlag = False
        self.nodeCounter = 0
        self.cycleFlag = False
        self.memoryAccessBytes = 0
        self.memoryWastedCycles = 0
        self.node1LoadedTiles = set([])
        self.node2LoadedTiles = set([])
        self.T = 32

    def setNext(self,PEArray,memory) -> None:
        self.PEArray = PEArray
        self.memory = memory
        
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
        # load in the PE metadata!: (post intersection fill) This doesn't happen at DRAM because there is no data storage at that level.
        coordCounter = 0
        if s1 == None:
            return 
        for x in s1.children:
            coordCounter += 1
            coordCounter += len(s1.children[x].children)
                
        for x in s2.children:
            coordCounter += 1
            coordCounter += len(s2.children[x].children)
                
        self.MemoryAccess(coordCounter * 2) # two bytes per coordinate
                
    
    
    def MemoryAccess(self, numBytes):
        self.memoryAccessBytes += numBytes
        self.memory.requestMemory(numBytes,"LLBI")
        self.memoryFlag = True
        
    def cycle(self) -> None:
        if self.endFlag or not self.inputFlag:
            return
    
        if self.memoryFlag:
            self.memoryWastedCycles += 1
            if self.memory.requestQueue["LLBI"]:
                return
            else:
                self.memoryFlag = False
        
        # If one of the streams runs out, we reset node 2 (B) and iterate to the next child of node 1 (A) (output A stationary)
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
                self.node1LoadedTiles = set([])
                self.node2LoadedTiles = set([])
                self.nodeCounter = 0
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                return
            self.streamOne = sorted(list(self.nodeOne.children[self.nodeOneIndices[self.nodeCounter]].children.keys()), reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children.keys()),reverse=True)
            self.nodeCounter += 1
        if self.skipto:
            # The skipTo technology essentially allows one stream to "catch up" to the other if it is multiple steps behind
            # It only allows for a given number, T, of skips. T represents the size of our comparator array. 
            counter = 0
            if self.streamOne[-1] < self.streamTwo[-1]:
                while (self.streamOne and self.streamTwo) and self.streamOne[-1] < self.streamTwo[-1] and counter < self.T:
                    self.streamOne.pop()
                    counter += 1
            elif self.streamOne[-1] > self.streamTwo[-1]:
                while (self.streamOne and self.streamTwo) and self.streamOne[-1] > self.streamTwo[-1] and counter < self.T:
                    self.streamTwo.pop()
                    counter += 1

            # if EOF or not equal, just go to the next cycle to load the next values.
            if not (self.streamOne and self.streamTwo) or self.streamOne[-1] != self.streamTwo[-1]:
                return
            
            k = self.streamOne.pop()
            node1Row = self.nodeOneIndices[self.nodeCounter-1]
            node = self.nodeOne.children[self.nodeOneIndices[self.nodeCounter-1]].children[k]
            node2List = self.nodeTwo.children[self.streamTwo.pop()].children.items()
            for node2 in node2List:
                # Since this coordinator is post intersection fill, this is the part where all of the PE data (metadata and scalars) are loaded into memory and sent because it is an
                # effectual calculation now. In a-stationary intersection, a PE tile can be used multiple times, so it is important to cache that data for later as well
                # NOTE: when a new LLB tile is loaded in, this is not the case. Even if the same tile is loaded in, the data memory is effectively wiped.
                if (node1Row,k) not in self.node1LoadedTiles:
                    self.node1LoadedTiles.add((node1Row,k))
                    byteCounter = 0
                    for x in node.children:
                        byteCounter += 2 # 2 bytes of coordinate data per I/J coordinate at the PE tile level (output stationary)
                        byteCounter += len(node.children[x].children) * 6 # 2 bytes for coordinate data per K coordinate at PE tile level, and 4 for its corresponding scalar int/float
                    self.MemoryAccess(byteCounter)
                # The reason for the swap here is that the columns of A and the Rows of B need to match to have an intersection, which is represented by K
                if (k,node2[0]) not in self.node2LoadedTiles:
                    self.node2LoadedTiles.add((k,node2[0]))
                    byteCounter = 0
                    for x in node2[1].children:
                        byteCounter += 2 # 2 bytes of coordinate data per I/J coordinate at the PE tile level (output stationary)
                        byteCounter += len(node2[1].children[x].children) * 6 # 2 bytes for coordinate data per K coordinate at PE tile level, and 4 for its corresponding scalar int/float
                    self.MemoryAccess(byteCounter)
                
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
        # S1 and S2 are PE tiles, meaning each individual intersection is given to a different PE!
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
        self.numEmptyCycles = 0
        self.memoryAccessBytes = 0
        self.T = 32

    
    def load(self, input: tuple) -> None:
        # put the offsets for the coordinates here. The final coords can be added to this to get the x and y of the output.
        # Note: We assume dataflow on-chip takes negligible time, and don't count it. 
        # In the paper, the full coordinates AS WELL AS the entire PE tile is sent through the NoC, so no data is pulled from memory.
        if input[0] != None:
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
        
    def MemoryAccess(self, numBytes):
        # This is actually unused, because the PE level intersectors get all the info they need from the LLB
        self.memoryAccessBytes += numBytes
    
    def cycle(self) -> None:   
        if self.endFlag or not self.inputFlag:
            self.numEmptyCycles += 1
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
                # Normally, memory access would be done here! Both nodes are fully loaded along with their scalars. Its easier
                # to count below, though so thats what I do.
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
            # Note: Memory access makes use of pre-intersection fills, meaning the data is loaded along with metadata.
            self.MemoryAccess((len(self.streamOne) * 2) * 3) # 2 bytes for coordinate, 4 bytes per coordinate for scalar values
            self.MemoryAccess((len(self.streamTwo) * 2) * 3) # 2 bytes for coordinate, 4 bytes per coordinate for scalar values
            return
        
        if self.skipto:
            counter = 0
            if self.streamOne[-1] < self.streamTwo[-1]:
                while (self.streamOne and self.streamTwo) and self.streamOne[-1] < self.streamTwo[-1] and counter < self.T:
                    self.streamOne.pop()
                    counter += 1
            elif self.streamOne[-1] > self.streamTwo[-1]:
                while (self.streamOne and self.streamTwo) and self.streamOne[-1] > self.streamTwo[-1] and counter < self.T:
                    self.streamTwo.pop()
                    counter += 1

            # if EOF or not equal, just go to the next cycle to load the next values.
            if not (self.streamOne and self.streamTwo) or self.streamOne[-1] != self.streamTwo[-1]:
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