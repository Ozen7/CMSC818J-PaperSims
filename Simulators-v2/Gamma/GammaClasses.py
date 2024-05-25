import numpy as np
from scipy.sparse import csr_matrix
from collections import deque
import threading

class Scheduler:

    def __init__(self,PEArray, Amatrix: csr_matrix, radix) -> None:
        self.PEArray = PEArray
        self.finishedJobs = set([])
        self.inProgressJobs = set([])
        self.jobs = deque([])
        self.Amatrix = Amatrix
        self.JobIDs = 0
        self.radix = radix
        self.jobLock = threading.Lock()
        
        for row in range(len(Amatrix.indptr)-1):
            self.jobs.extend(self.schedule(row,True,self.Amatrix.indptr[row],self.Amatrix.indptr[row+1]))
            
    def schedule(self,row,location,lptr,rptr):
        self.JobIDs += 1
        if rptr - lptr < self.radix:
            # location = output is final or partial, whether this input takes in A fibers and B fibers,
            # row indicates which row this belongs to, self.JobIDs is the id of the job, and the rest is just the data and indices
            if rptr > len(self.Amatrix.data):
                rptr = len(self.Amatrix.data)
            return (location, True, row, self.JobIDs, self.Amatrix.indices[lptr:rptr], self.Amatrix.data[lptr:rptr])
        else:
            # location = output is final or partial, this input takes in partial fibers and merges them
            # row indicates which row, self.jobIDs is the id of the job, the first list contains dependencies, the second is the 1.0 scalars
            retval = [(location, False, row, self.JobIDs, [], [1] * 64)]
            
            # this is the height of the tree. This is at least 2 (since if it was 1, it would be handled above)
            count = 0
            while self.radix**count <= rptr-lptr:
                count += 1
            
            # assign radix^(count-1) to each immediate subtree, when you run out, rebalance by taking radix^(count-2)
            # from full nodes and giving them to the remaining trees.
            distributionList = [0] * self.radix
            pointer = 0
            remainingValues = rptr-lptr
            for i in range(self.radix):
                if remainingValues >= self.radix**(count-1):
                    distributionList[i] = self.radix**(count-1)
                    remainingValues -= self.radix**(count-1)
                else:
                    # find a value in the distributionList that can "donate" merges to fill up the remaining slots
                    # there MUST be space, since we have more merges than one radix (the size of the tree is >= 2)
                    # if the remaining values is still large enough to assign more nodes, it will
                    while distributionList[pointer] <= (self.radix**(count-2) - remainingValues):
                        pointer += 1
                    
                    pointerloss = max(0,self.radix**(count-2) - remainingValues)
                    distributionList[pointer] -= pointerloss
                    remainingValues -= (self.radix**(count-2) - pointerloss)
                    distributionList[i] = self.radix**(count-2)
                
                # now, we use the distribution list in order to generate the subtrees
                iterator = 0
                for dist in distributionList:
                    o = self.schedule(row,False,iterator,iterator + dist) # summing up all the distances is equal to the number of fibers to be merged
                    iterator += dist
                    retval[0][4].append(o[3]) #add the job ID of this new job to the array of dependencies
                    retval += o #add the job to the ovrall return value
                
    def markCompleted(self, JobID):
        # add the job to the set of completed jobs - its already been saved to DRAM
        with self.jobLock:
            self.inProgressJobs.remove(JobID)
            self.finishedJobs.add(JobID)

    def reassign(self):
        # Now, we find the next job to assign to the PE that just finished working and is now switching to work on the task it accepted
        # before it finished (this task will be done after the one that was just started)
        # Note: the PE sends the fetch command to the FiberCache when it gets the output.
        # we use the lock in case two PEs try to get a new scheudle at the same time. Unlikely, but possible.
        with self.jobLock:
            for job in self.jobs:
                if job[3] in self.inProgressJobs:
                    continue
                if job[1]:
                    self.inProgressJobs.add(job[3])
                    return job
                else:
                    for x in job[4]:
                        if x not in self.finishedJobs:
                            continue
                    self.inProgressJobs.add(job[3])
                    return job
            else:
                # if we run through all jobs but can't find one that has all its requirements met, 
                # then this PE won't ever actually need to be used again. Remember - if a merge happens with partial matrices,
                # then the PE that did the last partial matrix merge can just pick up the job. This means that the PE wont ever have an 
                # opening to work again anyway, so we can just send None and tell the PE to stop computing
                return None
            
                
        
                
                
                
                
class FiberCache:
    # The core of the accelerator. Does the following:
    # Regular Caching, all memory that moves from or to memory and PEs gets cached in the FiberCache 
    # However, fetch and consume are added.
    # Fetch allows the fibercache to prefetch values into the fibercache, making it act as a buffer, and hiding overheads.
    # This basically makes sure that most read requests to the cache already have values preloaded in the cache, reducing cache misses
    # Consume removes a value that isn't going to be used from the Fibercache, allowing for memory to be more easily freed
    # Since the cache doesn't need to be backed by DRAM for the purposes of certain operations, we can just use the cache as a 
    # quick and dirty storage system to temporarily store stuff that the PE needs without having to send it back to DRAM
    # When a consume occurs, the written value inside the FiberCache is removed. When a write occurs that is not consumed, and the PE
    # changes the A row its working on, it is written back to memory. Read the paper closely to figure out what happens to data that
    # the PE has in the FiberCache when it is done operating on its row of A
    def __init__(self) -> None:
        pass


class PE:
    # Works along with the fibercache. Each PE works on a row of A at a time, merging partial matrices then mergeing them into complete
    # ones using a balanced, top-full tree scheduler
    
    def __init__(self, scheduler: "Scheduler", fiberCache: "FiberCache") -> None:
        self.scheduler = scheduler
        self.fiberCache = fiberCache
        self.endFlag = False
        
        self.nextFiberSet = self.scheduler.reassign()

        # nextFiberSet[1] is a boolean representing reading B fibers or partial fibers, and nextFiberSet[4] is a list of fibers to read
        if self.nextFiberSet[1]:
            self.fiberCache.fetch(self.nextFiberSet[4])

        
        self.currentFiberSet = None
        self.partialOrFinal = None
        self.outputPartialOrFinal = None
        self.row = None
        self.JobID = None
        self.inputFibers = None
        self.Adata = None
        
        
    def switchAndReassign(self):
    
        self.currentFiberSet = self.nextFiberSet
        if self.currentFiberSet != None:
            self.outputPartialOrFinal = self.currentFiberSet[0]
            self.inputPartialOrFinal = self.currentFiberSet[1]
            self.row = self.currentFiberSet[2]
            self.JobID = self.currentFiberSet[3]
            self.inputFibers = self.currentFiberSet[4]
            self.Adata = self.currentFiberSet[4]
        else:
            self.endFlag = True
            return
        
        if self.inputPartialOrFinal:
            for inp in range(len(self.inputFibers)):
                # setup B fibers data transfer from fiberCache into this input number
                # this fills in both the coords[inp] and the values[inp]
                self.fiberCache.read(self,inp,self.inputFibers[inp])
        else:
            for inp in range(len(self.inputFibers)):
                # setup partial fibers data transfer from fiberCache into this input number
                # this fills in both the coords[inp] and the values[inp]
                self.fiberCache.consume(self,inp,self.inputFibers[inp])
        
        
        self.nextFiberSet = self.scheduler.reassign()
        if self.nextFiberSet == None:
            return
        
        if self.nextFiberSet[1]:
            self.fiberCache.fetch(self.nextFiberSet[4])        
    
    def cycle(self):
        
        
        # there isn't anything mentioned in the paper, other than how many banks FiberCache has
        
        
        """
        1) If the output is running dry, stop compute
        2) have a set of booleans and queues just like SpArch where the boolean represents a queue running dry, and the queues are the fifos between layers (check if the register has a value or not, not if it is < 16)
        3) output coordinates are not merged, only one coordinate and one index is output every cycle. 
        4) We pop the top of the value fifo aligned with the index that is output this cycle, which is the value associated with the coordinate that is output this cycle
        5) We multiply the value with the scaling factor aligned with the index that is output
        6) We check whether the coordinate output this cycle is the same as last cycle. If it is, we sum them up and retain the combined value. 
            If not, the previously accumulated value and previous coordinate are output by the PE, either to be saved to memory (if final) or saved in fiberCache (if partial)
            Their places are taken by the coordinate and value that were output this cycle
        7) continue to do this until there is no more input left. output the final output coordinate and value, and move on to the next set of fibers to merge 
            (different coords, different values, different scaling factors). Contact the scheduler to let it know that the current job id is finished, and get the job after the new one assigned
        8) continue to do 7 until the scheduler runs out of jobs to give.    
        """
        
        
        if self.nextFiberSet == None:
            self.nextFiberSet = self.scheduler.reassign()

        
        
    
