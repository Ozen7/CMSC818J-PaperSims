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
from threading import Event

if __name__ == "__main__":
    LLB_TILE_SIZE = 10
    PE_TILE_SIZE = 5
    
    endFlag = True
    dramIntersector = DRAMIntersector(True)
    llbIntersector = LLBIntersector(True)
    peArray = PEArray(12)
    peIntersectorList = []
    for x in range(0,12):
        peIntersectorList.append(PEIntersector(True, LLB_TILE_SIZE, PE_TILE_SIZE,x))
    dramIntersector.setNext(llbIntersector)
    llbIntersector.setNext(peArray)
    peArray.setNext(peIntersectorList)
    
    
    
    gen = np.random.default_rng()
    data1 = gen.integers(1,100,10000)
    row1 = gen.integers(0,10000,10000)
    col1 = gen.integers(0,10000,10000)

    data2 = gen.integers(1,100,10000)
    row2 = gen.integers(0,10000,10000)
    col2 = gen.integers(0,10000,10000)
    i1 = coo_matrix((data1, (row1, col1)), shape=(10000, 10000))
    i2 = coo_matrix((data2, (row2, col2)), shape=(10000, 10000))

    input1 = coo_to_csf(i1,LLB_TILE_SIZE,LLB_TILE_SIZE,PE_TILE_SIZE,PE_TILE_SIZE,False)
    input2 = coo_to_csf(i2,LLB_TILE_SIZE,LLB_TILE_SIZE,PE_TILE_SIZE,PE_TILE_SIZE,True)

    #print("i1")
    #print_csf_tree(input1)
    #print("i2")
    #print_csf_tree(input2)

    #print(i1.toarray())
    #print(i2.toarray())
    dramIntersector.input(input1, input2)
    
    dIEvent = Event()
    dIEvent.set()
    Thread(target=dramIntersector.running,args=[dIEvent]).start()
    
    lIEvent = Event()
    lIEvent.set()
    Thread(target=llbIntersector.running,args=[lIEvent]).start()
    
    peAEvent = Event()
    peAEvent.set()
    Thread(target=peArray.running,args=[peAEvent]).start()
    
    peThreadList = [None] * 12
    peEventList = [None] * 12
    for x in range(0,12):
        peEventList[x] = Event()
        peEventList[x].set()
        Thread(target=peIntersectorList[x].running,args=[peEventList[x]]).start()
    count = 0
    while endFlag:
        count += 1
        dIEvent.clear()
        lIEvent.clear()
        peAEvent.clear()
        for event in peEventList:
            event.clear()

        if not dramIntersector.endFlag:
            dIEvent.wait()
        if not llbIntersector.endFlag:
            lIEvent.wait()
        if not peArray.endFlag:
            peAEvent.wait()
        for pe in range(0,12):
            if not peIntersectorList[pe].endFlag:
                peEventList[pe].wait()

        endFlag = not (dramIntersector.endFlag and llbIntersector.endFlag and peArray.endFlag)
        for x in range(0,12):
            endFlag = endFlag or not peIntersectorList[x].endFlag

    #r = []
    #c = []
    #v = []
    #for i,output in enumerate(peIntersectorList):
    #    for o in output.output:
    #        v.append(o[0])
    #        r.append(o[1])
    #        c.append(o[2])
    
    #est = coo_matrix((v,(r,c)),(10000,10000)).toarray()
    #actual = np.matmul(i1.toarray(),i2.toarray())
    #print(est)
    #print(actual)
    #print(np.equal(actual, est))
    #print(np.allclose(actual,est,0.0001,0.0001))
    print(count)
        
        
    
    
            
        