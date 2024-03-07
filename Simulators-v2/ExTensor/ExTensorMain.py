from ExTensorClasses import DRAMIntersector
from ExTensorClasses import LLBIntersector
from ExTensorClasses import PEArray
from ExTensorClasses import PEIntersector
from ExTensorClasses import CSFNode
from ExTensorClasses import coo_to_csf
from ExTensorClasses import print_csf_tree

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
        
        

    
    
            
        