from ExTensorClasses import DRAMIntersector
from ExTensorClasses import LLBIntersector
from ExTensorClasses import PEArray
from ExTensorClasses import PEIntersector

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

if __name__ == "__main__":
    endFlag = True
    dramIntersector = DRAMIntersector(40)
    llbIntersector = LLBIntersector()
    peArray = PEArray()
    peIntersectorList = []
    for x in range(0,12):
        peIntersectorList.append(PEIntersector())
    dramIntersector.setNext(llbIntersector)
    llbIntersector.setNext(peArray)
    peArray.setNext(peIntersectorList)

    
    while endFlag:
        endFlag = False
        endFlag = dramIntersector.cycle() or endFlag
        print(endFlag)
        endFlag = llbIntersector.cycle() or endFlag
        llbIntersector.load()
        print(endFlag)
        endFlag = peArray.cycle() or endFlag
        peArray.load()
        print(endFlag)
        for x in peIntersectorList:
            endFlag = x.cycle() or endFlag
            x.load()
            print(endFlag)
        
        

    
    
            
        