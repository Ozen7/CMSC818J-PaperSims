from ExTensorClasses import DRAMIntersector
from ExTensorClasses import LLBIntersector
from ExTensorClasses import PEArray
from ExTensorClasses import PEIntersector
from ExTensorClasses import CSFNode
from ExTensorClasses import coo_to_csf
from ExTensorClasses import print_csf_tree
import numpy as np
from scipy.sparse import coo_matrix
from threading import Thread

if __name__ == "__main__":
    endFlag = True
    dramIntersector = DRAMIntersector(True)
    llbIntersector = LLBIntersector(True)
    peArray = PEArray()
    peIntersectorList = []
    for x in range(0,12):
        peIntersectorList.append(PEIntersector())
    dramIntersector.setNext(llbIntersector)
    llbIntersector.setNext(peArray)
    peArray.setNext(peIntersectorList)
    
    gen = np.random.default_rng()
    data1 = gen.integers(1,100,3)
    row1 = gen.integers(0,4,3)
    col1 = gen.integers(0,4,3)

    data2 = gen.integers(1,100,3)
    row2 = gen.integers(0,4,3)
    col2 = gen.integers(0,4,3)
    i1 = coo_matrix((data1, (row1, col1)), shape=(4, 4))
    i2 = coo_matrix((data2, (row2, col2)), shape=(4, 4))

    input1 = coo_to_csf(i1,2,2,1,1,False)
    input2 = coo_to_csf(i2,2,2,1,1,True)

    print("i1")
    print_csf_tree(input1)
    print("i2")
    print_csf_tree(input2)

    print(i1.toarray())
    print(i2.toarray())
    dramIntersector.input(input1, input2)
    

    while endFlag:
        print("Iter")
        dI = Thread(target=dramIntersector.cycle)
        lI = Thread(target=llbIntersector.cycle)
        dI.start()
        lI.start()
        dI.join()
        lI.join()
        endFlag = not (dramIntersector.endFlag and llbIntersector.endFlag)

        '''
        print(endFlag)
        endFlag = not llbIntersector.endFlag or endFlag
        llbIntersector.load()
        print(endFlag)
        endFlag = not peArray.endFlag or endFlag
        peArray.load()
        print(endFlag)
        for x in peIntersectorList:
            endFlag = not x.endFlag or endFlag
            x.load()
            print(endFlag)
        '''
    for nodeset in peArray.inputBuffer:
        if nodeset[0]:
            print("NEWSET", "i1: ", nodeset[0].value, "i2: ", nodeset[1].value, "LLBRow: ", nodeset[2], "LLBCol: ", nodeset[3], "PERow: ", nodeset[4], "PECol: ", nodeset[5])
    print(np.matmul(i1.toarray(),i2.toarray()))
        
        

    
    
            
        