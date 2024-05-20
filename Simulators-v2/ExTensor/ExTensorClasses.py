from collections import deque
import time
import math
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
    def __init__(self, mode) -> None:
        self.mode = mode
        if self.mode != "NoMerge":
            self.streamOne = []
            self.streamTwo = []
        else:
            self.k = 0
        self.nodeOne = None
        self.nodeTwo = None
        self.j = 0
        self.LLBIntersector = None
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

        
    
    def setNext(self, llbi: 'LLBIntersector') -> None:
        self.LLBIntersector = llbi
    
    def input(self, s1: CSFNode, s2: CSFNode, matrix1Cols, LLBColSize):
        self.nodeOne = s1
        self.nodeTwo = s2
        self.j = 0
        self.nodeTwoIndices = sorted(list(s2.children.keys()))
        self.kMax = math.ceil(matrix1Cols/LLBColSize)
        # load in the LLB metadata!:
        coordCounter = 0
        for x in s1.children:
            coordCounter += 1
            for y in s1.children[x].children:
                coordCounter += 1
                
        for x in s2.children:
            coordCounter += 1
            for y in s2.children[x].children:
                coordCounter += 1
                
        self.MemoryAccess(coordCounter * 2) # two bytes per coordinate!
        # Normally, memory access would be done here! Both nodes are loaded up to the LLB metadata (due to post-intersection fill)
        # Its easier to count when individual rows are 
        
    def MemoryAccess(self, numBytes):
        self.memoryAccessBytes += numBytes
        self.memLatency = 70 + 15 + math.ceil(numBytes/128) #70 ns for fetching row, 15 for column latency, and then 128 bytes/cycle from then on, given 128gb/s and 1GhZ clock speed
        self.memoryFlag = True
        
        

    def cycle(self) -> None:
        if self.endFlag:                
            return
        
        # if there are things being pulled from memory, don't cycle.
        if self.memoryFlag:
            self.memoryWastedCycles += 1
            self.memLatency -= 1
            if self.memLatency <= 0:
                self.memoryFlag = False
            return

        # If we hit the limits of the intersection
        if self.mode != "NoMerge" and not (self.streamOne and self.streamTwo):
            if self.j >= len(self.nodeTwo.children):
                self.LLBIntersector.load(None, None, -1, -1)
                self.endFlag = True
                return
            
            self.streamOne = sorted(list(self.nodeOne.children.keys()),reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children[self.nodeTwoIndices[self.j]].children.keys()), reverse=True)
            self.j += 1
        elif self.mode == "NoMerge" and self.k > self.kMax:
            self.j += 1
            self.k = 0
            if self.j >= len(self.nodeTwo.children):
                self.LLBIntersector.load(None, None, -1, -1)
                self.endFlag = True
                return

        
        # intersection
        if self.mode == "Skip":
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
            
            # load node One's children one level down, and node Two's child two levels down. This is because both sides are in row major form
            node2Col = self.nodeTwoIndices[self.j-1]
            node2 = self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children[self.streamTwo.pop()]
            node1List = self.nodeOne.children[self.streamOne.pop()].children.items()
            for i in node1List:
                self.LLBIntersector.load(i[1], node2, i[0], node2Col )             
        elif self.mode == "NoMerge":
            if self.k in self.nodeOne.children and self.k in self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children:
                for i in self.nodeOne.children[self.k].children.items():
                    self.LLBIntersector.load(i[1], self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children[self.k], i[0], self.nodeTwoIndices[self.j - 1])
            self.k += 1         
        else:
            #Similar logic to the skipTo side above, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                node2Col = self.nodeTwoIndices[self.j-1]
                node2 = self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children[self.streamTwo.pop()]
                node1List = self.nodeOne.children[self.streamOne.pop()].children.items()
                for node in node1List:
                    self.LLBIntersector.load(node[1], node2, node[0], node2Col )          
            elif self.streamOne[-1] > self.streamTwo[-1]:
                self.streamTwo.pop()
            else:
                self.streamOne.pop()
    

class LLBIntersector:
    def __init__(self, mode, matrix1LLBColSize, PEColSize) -> None:
        self.mode = mode
        if self.mode != "NoMerge":
            self.streamOne = []
            self.streamTwo = []
        else:
            self.k = 0
            
        self.kMax = math.ceil(matrix1LLBColSize/ PEColSize)

        self.nodeOne = None
        self.nodeTwo = None
        self.PEArray = None
        self.endFlag = False
        self.nodeOneIndices = []
        self.inputBuffer = deque([])
        self.inputCoordBuffer = deque([])
        self.inputFlag = False
        self.memoryFlag = False
        self.i = 0
        self.memoryAccessBytes = 0
        self.memoryWastedCycles = 0
        self.node1LoadedTiles = set([])
        self.node2LoadedTiles = set([])
        self.T = 32

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
        self.memLatency = 5 + math.ceil(numBytes/128) # Column latency of 5 cycles + 128 bytes/cycle, given 128gb/s and 1GhZ clock speed and the required values already being in row buffer - latency for this is calculated once
        self.memoryFlag = True
        
    def cycle(self) -> None:
        if self.endFlag or not self.inputFlag:
            return
    
        if self.memoryFlag:
            self.memoryWastedCycles += 1
            self.memLatency -= 1
            if self.memLatency <= 0:
                self.memoryFlag = False
            return
        
        # If one of the streams runs out, we reset node 2 (B) and iterate to the next child of node 1 (A) (output A stationary)
        if self.mode != "NoMerge" and not (self.streamOne and self.streamTwo):
            if not (self.nodeOne or self.nodeTwo) or self.i >= len(self.nodeOne.children):
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
                self.memLatency = 70 # this is because of the time taken to bring row into row buffer
                self.i = 0
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                return
            self.streamOne = sorted(list(self.nodeOne.children[self.nodeOneIndices[self.i]].children.keys()), reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children.keys()),reverse=True)
            self.i += 1
        if (self.mode == "NoMerge" and self.k > self.kMax) or not (self.nodeOne or self.nodeTwo):
            self.k = 0
            self.i += 1
            if not (self.nodeOne or self.nodeTwo) or self.i >= len(self.nodeOne.children):
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
                self.memLatency = 70 # this is because of the time taken to bring row into row buffer
                self.i = 0
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                return


        if self.mode == "Skip":
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
            node1Row = self.nodeOneIndices[self.i-1]
            node = self.nodeOne.children[self.nodeOneIndices[self.i-1]].children[k]
            node2List = self.nodeTwo.children[self.streamTwo.pop()].children.items()
            for node2 in node2List:
                byteCounter = 0
                if (node1Row,k) not in self.node1LoadedTiles:
                    self.node1LoadedTiles.add((node1Row,k))
                    for x in node.children:
                        byteCounter += 2 # 2 bytes of coordinate data per I/J coordinate at the PE tile level (output stationary)
                        byteCounter += len(node.children[x].children) * 6 # 2 bytes for coordinate data per K coordinate at PE tile level, and 4 for its corresponding scalar int/float
                # The reason for the swap here is that the columns of A and the Rows of B need to match to have an intersection, which is represented by K
                if (k,node2[0]) not in self.node2LoadedTiles:
                    self.node2LoadedTiles.add((k,node2[0]))
                    for x in node2[1].children:
                        byteCounter += 2 # 2 bytes of coordinate data per I/J coordinate at the PE tile level (output stationary)
                        byteCounter += len(node2[1].children[x].children) * 6 # 2 bytes for coordinate data per K coordinate at PE tile level, and 4 for its corresponding scalar int/float
                if byteCounter:
                    self.MemoryAccess(byteCounter)
                
                self.PEArray.load(node, node2[1], self.row, self.col, node1Row, node2[0])
        elif self.mode == "NoMerge":
            byteCounter = 0
            if self.k in self.nodeTwo.children and self.k in self.nodeOne.children[self.nodeOneIndices[self.i]].children:
                for j in self.nodeTwo.children[self.k].children.items():
                    self.PEArray.load(self.nodeOne.children[self.nodeOneIndices[self.i]].children[self.k], j[1], self.row, self.col, self.nodeOneIndices[self.i], j[0])
                    for x in self.nodeTwo.children[self.k].children[j[0]].children:
                        byteCounter += 2 # 2 bytes of coordinate data per K coordinate at the PE tile level
                        len(self.nodeTwo.children[self.k].children[j[0]].children[x].children) * 6 # 2 bytes of data for each J coordinate, 4 more for scalar data
                    
                    for x in self.nodeOne.children[self.nodeOneIndices[self.i]].children[self.k].children:
                        byteCounter += 2
                        len(self.nodeOne.children[self.nodeOneIndices[self.i]].children[self.k].children[x].children) * 6
            if byteCounter:
                self.MemoryAccess(byteCounter)     
            self.k += 1
        else:
            #Similar logic to skipTo, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                k = self.streamOne.pop()
                node1Row = self.nodeOneIndices[self.i-1]
                node = self.nodeOne.children[self.nodeOneIndices[self.i-1]].children[k]
                node2List = self.nodeTwo.children[self.streamTwo.pop()].children.items()
                for node2 in node2List:
                    byteCounter = 0
                    if (node1Row,k) not in self.node1LoadedTiles:
                        self.node1LoadedTiles.add((node1Row,k))
                        for x in node.children:
                            byteCounter += 2 # 2 bytes of coordinate data per I/J coordinate at the PE tile level (output stationary)
                            byteCounter += len(node.children[x].children) * 6 # 2 bytes for coordinate data per K coordinate at PE tile level, and 4 for its corresponding scalar int/float
                    # The reason for the swap here is that the columns of A and the Rows of B need to match to have an intersection, which is represented by K
                    if (k,node2[0]) not in self.node2LoadedTiles:
                        self.node2LoadedTiles.add((k,node2[0]))
                        for x in node2[1].children:
                            byteCounter += 2 # 2 bytes of coordinate data per I/J coordinate at the PE tile level (output stationary)
                            byteCounter += len(node2[1].children[x].children) * 6 # 2 bytes for coordinate data per K coordinate at PE tile level, and 4 for its corresponding scalar int/float
                    if byteCounter:
                        self.MemoryAccess(byteCounter)
                    
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
    def __init__(self, mode, llbTileSize, PETileSize, num) -> None:
        self.mode = mode
        if self.mode != "NoMerge":
            self.streamOne = []
            self.streamTwo = []
        else:
            self.k = -1
            
        self.kMax = PETileSize
        self.nodeOne = None
        self.nodeTwo = None
        self.output = []
        self.mode = mode
        self.endFlag = False
        self.inputFlag = False
        self.nodeOneIndices = []
        self.i = 0
        self.j = 0
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
    
    def cycle(self) -> None: 
        if self.endFlag or not self.inputFlag:
            self.numEmptyCycles += 1
            return
        
        # If one of the streams runs out, we reset node 1 (A) and iterate to the next child of node 2 (B) (output B stationary)
        if self.mode != "NoMerge" and not (self.streamOne and self.streamTwo):
            if self.nodeOne and self.nodeTwo and self.j >= len(self.nodeTwo.children) and self.i < len(self.nodeOne.children):
                self.i += 1
                self.j = 0
            elif not (self.nodeOne or self.nodeTwo) or self.j >= len(self.nodeTwo.children):
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
                self.j = 0
                self.i = 1
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                self.nodeTwoIndices = sorted(list(self.nodeTwo.children.keys()))
                return
            self.j += 1
            self.streamOne = sorted(list(self.nodeOne.children[self.nodeOneIndices[self.i-1]].children.keys()),reverse=True)
            self.streamTwo = sorted(list(self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children.keys()),reverse=True)
            # Note: Memory access makes use of pre-intersection fills, but all data is sent from the LLB intersector, so no DRAM access is necessary
            return
        elif self.mode == "NoMerge" and self.k > self.kMax or not (self.nodeOne or self.nodeTwo) or self.k == -1:
            if self.nodeOne and self.nodeTwo and self.j >= len(self.nodeTwo.children) and self.i < len(self.nodeOne.children):
                self.i += 1
                self.j = 0
            elif not (self.nodeOne or self.nodeTwo) or self.j >= len(self.nodeTwo.children):
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
                self.j = 0
                self.i = 1
                self.nodeOneIndices = sorted(list(self.nodeOne.children.keys()))
                self.nodeTwoIndices = sorted(list(self.nodeTwo.children.keys()))
                return
            self.k = 0
            self.j += 1
            
        
        if self.mode == "Skip":
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

            node1Row = self.nodeOneIndices[self.i-1]
            node2Col = self.nodeTwoIndices[self.j-1]
            node1 = self.nodeOne.children[self.nodeOneIndices[self.i-1]].children[self.streamOne.pop()]
            node2 = self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children[self.streamTwo.pop()]
            self.output.append((node1.value * node2.value, self.row + node1Row, self.col + node2Col))
        elif self.mode == "NoMerge":
            #print("Iter " + str(self.id) + " " + str(self.i) + " " + str(self.j) + " " + str(self.k))
            if self.k in self.nodeOne.children[self.nodeOneIndices[self.i-1]].children and self.k in self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children:
                n1 = self.nodeOne.children[self.nodeOneIndices[self.i-1]].children[self.k]
                n2 = self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children[self.k]
                #print("Output " + str(self.id) + " " + str((n1.value * n2.value, self.row + self.nodeOneIndices[self.i-1],self.col + self.nodeTwoIndices[self.j-1])))
                self.output.append((n1.value * n2.value, self.row + self.nodeOneIndices[self.i-1],self.col + self.nodeTwoIndices[self.j-1]))
            self.k += 1
        else:
            #Similar logic to the skipTo side above, but instead of skipping to the correct value we go through one by one.
            if self.streamOne[-1] == self.streamTwo[-1]:
                node1Row = self.nodeOneIndices[self.i-1]
                node2Col = self.nodeTwoIndices[self.j-1]
                node1 = self.nodeOne.children[self.nodeOneIndices[self.i-1]].children[self.streamOne.pop()]
                node2 = self.nodeTwo.children[self.nodeTwoIndices[self.j-1]].children[self.streamTwo.pop()]
                self.output.append((node1.value * node2.value, self.row + node1Row, self.col + node2Col))
            elif self.streamOne[-1] > self.streamTwo[-1]:
                self.streamTwo.pop()
            else:
                self.streamOne.pop()
                
                
# This is extra code made specifically to test ExTensor without hierarchical intersection
# Basically, instead of finding matches that actually amount to something in the DRAM and LLB layers, we just 
# iterate through coordinates of PE tiles and assign them to individual PEs

# Step 1) Rewrite CSF code so it only has PE level tiles:

def coo_to_csf_small(coo_matrix, llbxsize, llbysize, pexsize, peysize, swap):
    root = CSFNode(None)

    for i, j, data in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
        current_node = root
        
        
        c = i // pexsize
        d = j // peysize
        
        
        e = (i - (c * pexsize))
        f = (j - (d * peysize))
        
        if swap:
            temp = e
            e = f
            f = temp

        # Traverse the tree based on the coordinates
        # NOTICE: b comes before a because THE LLB TILES ARE COLUMN MAJOR, NOT ROW MAJOR. BECAUSE THE MULTIPLICATION IS B-STATIONARY, THE MATRIX MUST BE COL MAJOR FOR THE FIRST TWO COORDINATES
        # WHEN DOING A-STATIONARY, WE WANT BOTH AS ROW-MAJOR, AND WHEN OUTPUT STATIONARY WE WANT A AS ROW MAJOR, AND B AS COL MAJOR
        for coord in [c, d, e, f]:
            if coord not in current_node.children:
                current_node.children[coord] = CSFNode(None)
            current_node = current_node.children[coord]

        # Set the leaf node value to the data value
        if current_node.value != None:
            current_node.value += data
        else:
            current_node.value = data

    return root




# Step 2) Write a new "intersector" that acts as a memory loader
# Basically, it iterates through the tree in a "A-stationary" manner. It loads in a value from row A, and then loads in the corresponding
# Row from B. It then sends these values into the PE array, where they get distributed to PEs that work the same way.
# if there is a zero-valued side, we don't need to load anything, but it wastes a cycle.
# More or less, we check every possibility one by one instead of skipping all 0-valued PE tiles.
class DRAMBaselineIntersector:
    def __init__(self) -> None:
        pass


class LLBBaselineIntersector:
    def __init__(self, PeRow, PeCol) -> None:
        self.i = 0
        self.k = 0
        self.j = 0
        
        self.PeCol = PeCol # the number of columns in a PE tile
        self.PeRow = PeRow # the number of rows in a PE tile
        
        self.nodeOne = None
        self.nodeTwo = None
        
        self.PEArray = None
        
        self.endFlag = False
        self.inputFlag = False
        self.memoryFlag = False
        
        self.nodeOneIndices = []
        self.inputBuffer = deque([])
        self.inputCoordBuffer = deque([])

        self.memoryAccessBytes = 0
        self.memoryWastedCycles = 0
        
    def setNext(self,PEArray) -> None:
        self.PEArray = PEArray
        
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def load(self, s1: CSFNode|None, s2: CSFNode|None, Matrix1Row, Matrix2Row, Matrix2Col) -> None:
        self.nodeOne = s1
        self.nodeTwo = s2
        
        self.nodeOneKeys = s1.keys
        
        self.iMax = math.ceil(Matrix1Row/self.PeRow) # total number of PE tile rows
        self.kMax = math.ceil(Matrix2Row/self.PeRow)
        self.jMax = math.ceil(Matrix2Col/self.PeCol) # total number of PE tile columns

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
        self.memLatency = 5 + math.ceil(numBytes/128) # Column latency of 5 cycles + 128 bytes/cycle, given 128gb/s and 1GhZ clock speed and the required values already being in row buffer - latency for this is calculated once
        self.memoryFlag = True
        
    def cycle(self) -> None:
        if self.endFlag or not self.inputFlag:
            return
    
        if self.memoryFlag:
            self.memoryWastedCycles += 1
            self.memLatency -= 1
            if self.memLatency <= 0:
                self.memoryFlag = False
            return
        
        
        if self.k > self.kMax:
            self.i += 1
            self.k = 0
        
        if self.i > self.iMax:
            self.endFlag = True
            return
        
        
        # the critical difference here is that the data does not iterate by nonzero values, and goes 1 by 1
        # Basically, instead of skipping to the next valid K coordinate in either stream and trying to find matches, we just 
        # iterate K one by one each cycle and match any nonzeroes
        
        # Basically, check if the row k exists, and whether the point i,k exists. If they do, match all nonzeros. If not, iterate k.
        # In the base dataflow, we match k dynamically by going through all nonzero k's and finding matches. Here, we just iterate k.
        byteCounter = 0
        if self.k in self.nodeTwo.children and self.i in self.nodeOne.children and self.k in self.nodeOne.children[self.i]:
            for j in self.nodeTwo.children[self.k]:
                self.PEArray.load(self.nodeOne.children[self.i].children[self.k], self.nodeTwo.children[self.k].children[self.j], 0, 0, self.i, self.j)
                for x in self.node2.children[self.k].children[self.j].children:
                    byteCounter += 2 # 2 bytes of coordinate data per K coordinate at the PE tile level
                    len(self.node2.children[self.k].children[self.j].children[x]) * 6 # 2 bytes of data for each J coordinate, 4 more for scalar data
                
                for x in self.node1.children[self.i].children[self.k].children:
                    byteCounter += 2
                    len(self.node1.children[self.i].children[self.k].children[x]) * 6
        if byteCounter:
            self.MemoryAccess(byteCounter)
                  
        self.k += 1
        
        


# Step 3) Create new PE that also works by iterating K and not only looking at nonzero values.

class PEBaselineIntersector:
    def __init__(self, PeRow, PeCol) -> None:
        self.i = 0
        self.k = 0
        self.j = 0
        
        self.PeCol = PeCol # the number of columns in a PE tile
        self.PeRow = PeRow # the number of rows in a PE tile
        
        self.nodeOne = None
        self.nodeTwo = None
        
        self.PEArray = None
        
        self.endFlag = False
        self.inputFlag = False
        self.memoryFlag = False
        
        self.nodeOneIndices = []
        self.inputBuffer = deque([])
        self.inputCoordBuffer = deque([])

        self.memoryAccessBytes = 0
        self.memoryWastedCycles = 0
        

    def setNext(self,PEArray) -> None:
        self.PEArray = PEArray
        
    def running(self, event):
        while not self.endFlag:
            time.sleep(0.0001)  
            if not event.isSet():
                self.cycle()
                event.set()
        event.set()
    
    def load(self, s1: CSFNode|None, s2: CSFNode|None, Matrix1Row, Matrix2Row, Matrix2Col) -> None:
        self.nodeOne = s1
        self.nodeTwo = s2
        
        self.iMax = math.ceil(Matrix1Row/self.PeRow) # total number of PE tile rows
        self.kMax = math.ceil(Matrix2Row/self.PeRow)
        self.jMax = math.ceil(Matrix2Col/self.PeCol) # total number of PE tile columns

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
        self.memLatency = 5 + math.ceil(numBytes/128) # Column latency of 5 cycles + 128 bytes/cycle, given 128gb/s and 1GhZ clock speed and the required values already being in row buffer - latency for this is calculated once
        self.memoryFlag = True
        
    def cycle(self) -> None:
        if self.endFlag or not self.inputFlag:
            return
    
        if self.memoryFlag:
            self.memoryWastedCycles += 1
            self.memLatency -= 1
            if self.memLatency <= 0:
                self.memoryFlag = False
            return
        
        
        if self.k > self.kMax:
            self.i += 1
            self.k = 0
        
        if self.i > self.iMax:
            self.endFlag = True
            return
        
        
        # the critical difference here is that the data does not iterate by nonzero values, and goes 1 by 1
        # Basically, instead of skipping to the next valid K coordinate in either stream and trying to find matches, we just 
        # iterate K one by one each cycle and match any nonzeroes
        
        # Basically, check if the row k exists, and whether the point i,k exists. If they do, match all nonzeros. If not, iterate k.
        # In the base dataflow, we match k dynamically by going through all nonzero k's and finding matches. Here, we just iterate k.
        byteCounter = 0
        if self.k in self.nodeTwo.children and self.i in self.nodeOne.children and self.k in self.nodeOne.children[self.i]:
            for j in self.nodeTwo.children[self.k]:
                self.PEArray.load(self.nodeOne.children[self.i].children[self.k], self.nodeTwo.children[self.k].children[self.j], 0, 0, self.i, self.j)
                for x in self.node2.children[self.k].children[self.j].children:
                    byteCounter += 2 # 2 bytes of coordinate data per K coordinate at the PE tile level
                    len(self.node2.children[self.k].children[self.j].children[x]) * 6 # 2 bytes of data for each J coordinate, 4 more for scalar data
                
                for x in self.node1.children[self.i].children[self.k].children:
                    byteCounter += 2
                    len(self.node1.children[self.i].children[self.k].children[x]) * 6
        if byteCounter:
            self.MemoryAccess(byteCounter)
                  
        self.k += 1