import numpy as np
from scipy.sparse import csr_matrix
from collections import deque
import threading

class Scheduler:

    def __init__(self,PEArray, Amatrix: csr_matrix, radix) -> None:
        self.PEArray = PEArray
        self.finishedJobs = set([])
        self.inProgressJobs = set([])
        self.jobs = []
        self.Amatrix = Amatrix
        self.JobIDs = 0
        self.radix = radix
        self.jobLock = threading.Lock()
        
        for row in range(len(Amatrix.indptr)-1):
            # "huge" assumption made here: scheduler only schedules rows with nonzeros in them
            if self.Amatrix.indptr[row] - self.Amatrix.indptr[row+1] == 0:
                continue
            self.jobs.extend(self.schedule(row,True,self.Amatrix.indptr[row],self.Amatrix.indptr[row+1]))
            #print(self.jobs)
            
    # row  = which row is being scheduled, location = whether the output is partial or not, lptr and rptr are left and right pointers
    def schedule(self,row,location,lptr,rptr):
        self.JobIDs += 1
        if rptr - lptr <= self.radix:
            # location = output is final or partial, whether this input takes in A fibers and B fibers,
            # row indicates which row this belongs to, self.JobIDs is the id of the job, and the rest is just the data and indices
            if rptr > len(self.Amatrix.data):
                rptr = len(self.Amatrix.data)
            return [(location, [True for _ in range(rptr-lptr)], row, self.JobIDs, self.Amatrix.indices[lptr:rptr], list(self.Amatrix.data[lptr:rptr]))]
        else:
            # location = output is final or partial, this input takes in partial fibers and merges them
            # row indicates which row, self.jobIDs is the id of the job, the first list contains dependencies, the second is the 1.0 scalars
            retval = [(location, [], row, self.JobIDs, [], [])]
            
            # this is the height of the tree. This is at least 2 (since if it was 1, it would be handled above)
            count = 1
            while self.radix**count < rptr-lptr:
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
                    while distributionList[pointer] < (self.radix**(count-1)-remainingValues):
                        pointer += 1
                    
                    # if the remaining values are enough to cover it, they do
                    pointerloss = max(0,self.radix**(count-1) - remainingValues)
                    # donate pointerloss values from distributionList
                    distributionList[pointer] -= pointerloss
                    # donate the remainder from remainingValues
                    remainingValues -= (self.radix**(count-1) - pointerloss)
                    distributionList[i] = self.radix**(count-1)
                    
            # ensures that there is a minimal number of merge operations at the lowest level
            if count == 2:
                merged_list = []
                current_sum = distributionList[0]
                for num in distributionList[1:]:
                    if current_sum + num <= self.radix:
                        current_sum += num
                    else:
                        merged_list.append(current_sum)
                        current_sum = num

                merged_list.append(current_sum)
                
                distributionList = merged_list
            
            # now, we use the distribution list in order to generate the subtrees
            iterator = lptr
            for dist in distributionList:
                # if dist is 0, we don't include the job
                if dist == 0:
                    continue
                # if the dist is only 1, we include it as an individual value to be merged
                elif dist == 1:
                    retval[0][1].append(True)
                    retval[0][4].append(int(self.Amatrix.indices[iterator]))
                    retval[0][5].append(int(self.Amatrix.data[iterator]))
                    iterator += dist
                    continue
                # otherwise, we recursively schedule
                o = self.schedule(row,False,iterator,iterator + dist) # summing up all the distances is equal to the number of fibers to be merged
                retval[0][1].append(False)
                retval[0][4].append(int(o[0][3])) # add the job ID of this new job to the array of dependencies
                retval[0][5].append(1)
                retval.extend(o) #add the job to the ovrall return value
                iterator += dist
            return retval

                
    def markCompleted(self, JobID):
        # add the job to the set of completed jobs - its already been saved to DRAM
        if JobID == None:
            return
        with self.jobLock:
            self.inProgressJobs.remove(JobID)
            self.finishedJobs.add(JobID)


    def reassign(self):
        # Now, we find the next job to assign to the PE that just finished working and is now switching to work on the task it accepted
        # before it finished (this task will be done after the one that was just started)
        # Note: the PE sends the fetch command to the FiberCache when it gets the output.
        # we use the lock in case two PEs try to get a new scheudle at the same time. Unlikely, but possible.
        #print(self.jobs)
        #print(self.finishedJobs)

        with self.jobLock:
            j = None
            for i,job in enumerate(self.jobs):
                t = True
                for x in job[1]:
                    t = t and x
                if t:
                    j = self.jobs.pop(i)
                    break
                else:
                    available = True
                    for t in range(len(job[4])):
                        if not job[1][t] and job[4][t] not in self.finishedJobs:
                            available = False
                    if available:
                        j = self.jobs.pop(i)
                        break
            
            # if we run through all jobs but can't find one that has all its requirements met, 
            # then this PE won't ever actually need to be used again. Remember - if a merge happens with partial matrices,
            # then the PE that did the last partial matrix merge can just pick up the job. This means that the PE wont ever have an 
            # opening to work again anyway, so we can just send None and tell the PE to stop computing
            if j != None:
                self.inProgressJobs.add(j[3])
            return j
            
                
        
           
         
# For greater realism, add a "job" part to each request, so when the fiberCache picks it back up it knows which job this 
# memory request is fulfilling.
class Memory:
    def __init__(self, numChannels, peakBandwidthPerChannel) -> None:
        self.numChannels = numChannels
        self.channelQueues = [deque() for _ in range(numChannels)]
        self.channelCurrent = [(None,0,0) for _ in range(numChannels)]
        self.BPC = peakBandwidthPerChannel
        self.TotalMemoryPulled = 0
        self.NumCyclesInUse = 0 #use this for memory efficiency when in use.
    
    def requestData(self,loader, ID, channel, amount):
        # numCycles includes the row and column latencies
        self.channelQueues[channel].append([loader,amount,ID]) # (loader, number of bytes)
    
    def cycle(self):
        inuse = False
        for i,channel in enumerate(self.channelCurrent):
            if channel[1] > 0:
                #print("Memory Push: ", channel)
                inuse = True
                sent = min(channel[1],self.BPC)
                self.TotalMemoryPulled += sent
                channel[0].inputQueue.append((channel[2],sent))
                channel[1] -= sent
            elif channel[1] <= 0 and self.channelQueues[i]:
                self.channelCurrent[i] = self.channelQueues[i].popleft()
        if inuse:
            self.NumCyclesInUse += 1              
                
  
# Note: FiberCache has to run after the PEs are done. This is because it could finish loading a fiber onto cache, and then 
# remove it before a PE has the chance to lock it for its own use!              
class FiberCache:
    def __init__(self, fiberCacheSize, B:csr_matrix) -> None:
        self.B = B
        self.remainingSize = fiberCacheSize
        self.fetchQueue = deque([])
        self.inputQueue = deque([])
        self.fetchSet = set([])
        self.fetchAmounts = {}
        self.BfiberSizes = {} # B fibers are stored by row value
        self.PartialValueSizes = {} # partial fibers are stored by Job ID that generated them.
        self.PartialValuesCoords = {}
        self.PartialValuesVals = {}
        self.PartialDRAM = {}
        self.BfiberPriority = {}
        self.Memory = None
        
        self.BFiberlockSet = set([])
        self.partialFiberLockSet = set([])
        self.outputMatrixIndptrs = []
        self.outputMatrixData = []
        self.outputMatrixCoords = []
        
        self.memoryUsage = 0
        self.memoryID = 0
    
    def setMemory(self,M:"Memory"):
        self.Memory = M
    
    def spillLines(self, amount):
        # first look through unlocked B fibers and replace lowest prio first
        # then, look at unlocked partial fibers and replace the lowest prio there
        # if everything is locked, we can't choose a replacement, and have to wait.
        possible = 0
        sorted_prio = [k for k,v in sorted(self.BfiberPriority.items(), key=lambda item: item[1],reverse=True)]
        while sorted_prio and amount > possible:
            row = sorted_prio.pop()
            # if the fiber is not locked or currently being fetched:
            if row not in self.BFiberlockSet and (True, row) not in self.fetchSet:
                possible += self.BfiberSizes[row]
        
        if amount > possible:
            # we start at the smallest partial values and go up. There is no priority for partial fibers
            sorted_size = [k for k,v in sorted(self.PartialValueSizes.items(), key=lambda item: item[1],reverse=True)]
            while sorted_size and amount > possible:
                JobID = sorted_size.pop()
                if JobID not in self.partialFiberLockSet and (False,JobID) not in self.fetchSet:
                    possible += self.PartialValueSizes[JobID]
        
        # if we are unable to muster up the right amount of space, we return false
        if amount > possible:
            return False
        
        sorted_prio = [k for k,v in sorted(self.BfiberPriority.items(), key=lambda item: item[1],reverse=True)]
        while sorted_prio and amount > 0:
            row = sorted_prio.pop()
            # if the fiber is not locked or currently being fetched:
            if row not in self.BFiberlockSet and (True, row) not in self.fetchSet:
                self.BfiberPriority.pop(row)
                s = self.BfiberSizes.pop(row)
                amount -= s
                self.remainingSize += s
        
        if amount > 0:
            # we start at the smallest partial values and go up. There is no priority for partial fibers
            sorted_size = [k for k,v in sorted(list(self.PartialValuesVals.items()), key=lambda item: item[1],reverse=True)]
            while sorted_size and amount > 0:
                JobID = sorted_size.pop()
                if JobID not in self.partialFiberLockSet and (False,JobID) not in self.fetchSet:
                    size = self.PartialValueSizes[JobID]
                    data = self.PartialValuesVals.pop(JobID)
                    coords = self.PartialValuesCoords.pop(JobID)
                    self.PartialDRAM[JobID] = (size,data,coords)
                    amount -= size
                    self.remainingSize += size
        return True
                       
    
    # operates to fetch fibers from dram, and sets priority
    def fetch(self, partialOrFinal, row):
        if partialOrFinal:
            if row in self.BfiberSizes or (True,row) in self.fetchSet:
                self.BfiberPriority[row] += 1
            else:
                # adds to the fetchQueue, where it will be scheduled to be pulled from memory.
                self.BfiberSizes[row] = (self.B.indptr[row+1] - self.B.indptr[row]) * 8 # 4 for data, 4 for coord
                self.BfiberPriority[row] = 1
                if self.BfiberSizes[row] != 0:
                    self.fetchQueue.append((True,row)) # True = B fiber, False = Partial Fiber
                    self.fetchSet.add((True,row))
        else:
            if row in self.PartialValuesVals or (False,row) in self.fetchSet:
                return
            else:
                (l,v,c) = self.PartialDRAM.pop(row)
                self.PartialValuesCoords[row] = c
                self.PartialValuesVals[row] = v
    
    def save(self, partialOrFinal, JobID, row, coords, vals):
        #print("SAVING", partialOrFinal, JobID, row, coords, vals)
        if partialOrFinal == None or (len(coords) == 0 and partialOrFinal):
            return
        if partialOrFinal:
            # saving to DRAM in CSR format, might add memory accesses here!
            self.outputMatrixIndptrs.append((row,len(coords)))
            self.outputMatrixCoords += coords
            self.outputMatrixData += vals
        else:
            l = len(coords) * 2 * 4 # data size of the coordinates and values, each number is 4 bytes
            
            self.PartialValueSizes[JobID] = l # saves the size in case we need to spill it or pull it

            # spill lines if necessary to free up space on cache. The dirty bit is implicit here. If it is spilled, 
            # the value is saved back to DRAM
            able = True
            if not self.remainingSize > l:
                able = self.spillLines(l - self.remainingSize)
            
            # if we can't spill enough lines to fit in the partial fiber, we save directly to DRAM
            if not able:
                self.PartialDRAM[JobID] = (l,vals,coords)
                return
            
            self.remainingSize -= l # row number is not needed to be stored, coord and value are each 4 bytes
            
            self.PartialValuesVals[JobID] = vals # saves the actual values so they can be fetched when needed by the PEs
            self.PartialValuesCoords[JobID] = coords
                
    
    # this ONLY runs on B fibers that have already been fetched and locked! 
    # it is great cause for error if this is not the case
    def read(self, PE: "PE", inputWayNumber, rowNumber):
        # load the values in this row into the PE way
        PE.coordFifos[inputWayNumber] = deque(self.B.indices[self.B.indptr[rowNumber]: self.B.indptr[rowNumber+1]])
        PE.valueFifos[inputWayNumber] = deque(self.B.data[self.B.indptr[rowNumber]:self.B.indptr[rowNumber+1]])
     
    
    # again, this only runs if all the partial fibers are in the cache, and have been locked already
    # we do NOT consider transfer throughput or latency between the FiberCache and the PE, its not something the 
    # paper really cares about, so we won't either
    def consume(self, PE: "PE", inputWayNumber, JobID):
        
        PE.coordFifos[inputWayNumber] = deque(self.PartialValuesCoords[JobID])
        PE.valueFifos[inputWayNumber] = deque(self.PartialValuesVals[JobID])
        
        # now, we delete these values from cache since they have no reuse
        self.PartialValuesCoords.pop(JobID)
        self.PartialValuesVals.pop(JobID)
        self.remainingSize += self.PartialValueSizes.pop(JobID)
    
    def lockFiber(self, partialOrFinal, fiber):
        if partialOrFinal:
            self.BFiberlockSet.add(fiber)
        else:
            self.partialFiberLockSet.add(fiber)
    
    def unlockFibers(self, partialOrFinal, fibers):
        for i,x in enumerate(partialOrFinal):
            if x:
                if fibers[i] in self.BFiberlockSet:
                    self.BFiberlockSet.remove(fibers[i])
            else:
                if fibers[i] in self.partialFiberLockSet:
                    self.partialFiberLockSet.remove(fibers[i])
    
    def checkIfFiberInCache(self, partialOrFinal, ID):
        if partialOrFinal:
            # if the value is in the cache and isn't being fetched
            if ID in self.BfiberSizes and (True,ID) not in self.fetchSet:
                return True
            if (True,ID) not in self.fetchSet:
                self.fetch(True,ID)
            return False
        else:
            # absolutely must be in partialValueSizes no matter what.
            if ID in self.PartialValuesVals and (False,ID) not in self.fetchSet:
                return True
            if (False,ID) not in self.fetchSet:
                self.fetch(False,ID)
            return False
        
    
    def cycle(self):
        '''
        Then, we iterate through the fetchQueue and send requests to memory. These can come in two forms:
        1) pulling B fibers (rows) from memory
        2) pulling partial fibers from memory
        
        First, we check over all incoming values and iterate through all of those,
        1) keep track of which request each value returned from DRAM associates with. 
        2) If we finish loading in a given value, we check it off on the fetchSet
            i) if a PE notices that the value has been checked off, it will load it individually and lock it
            i) it is locked because we fetched these values with the full intention of using them very very soon. 
               The only unlocked values in the Cache should be those that have already been used and aren't going to be used 
               either. 
        3) This should cause the checkIfFiberInCache method to return true for that given value.
        '''

        for _ in range(len(self.fetchQueue)):
            req = self.fetchQueue[0]
            if req[0]: # pulling B fibers
                size = self.BfiberSizes[req[1]]
            else:
                size = self.PartialValueSizes[req[1]]


            possible = True
            if size > self.remainingSize:
                possible = self.spillLines(size)
            
            if not possible:
                break # we need to wait in order to pull the next line

            self.memoryID += 1  
            self.fetchQueue.popleft()
            # pulling values from dram, split into 8 byte chunks to make full use of channel-level parallelism
            O = 0
            C = 0
            for x in range(0,size//8):
                ms = 8
                self.Memory.requestData(self,self.memoryID,C,ms)
                self.memoryUsage += ms
                # increment offset. If we move to a new channel, increment the channel number mod 8
                O += 8
                O %= 64
                if O < 8:
                    C += 1
                    C %= 8
            
            self.fetchAmounts[self.memoryID] = [req, size]
        # manage data pulled from DRAM this cycle
        for _ in range(len(self.inputQueue)):
            inp = self.inputQueue.popleft()
            self.fetchAmounts[inp[0]][1] -= inp[1]
            if self.fetchAmounts[inp[0]][1] <= 0:
                self.fetchSet.remove(self.fetchAmounts[inp[0]][0])
                self.fetchAmounts.pop(inp[0])


class PE:
    # Works along with the fibercache. Each PE works on a row of A at a time, merging partial matrices then mergeing them into complete
    # ones using a balanced, top-full tree scheduler
    
    def __init__(self, radix) -> None:
        self.endFlag = False
        self.dataNotInCacheFlag = False
        self.radix = radix
        self.total = 0 


        # this is NOT a merge tree! each cycle, the tree of compute units chooses a minimum value, and leaves the other queues alone
        # that is, no data is stored outside of the input fifos.
        self.coordFifos = [deque([]) for _ in range(self.radix)]
        self.valueFifos = [deque([]) for _ in range(self.radix)]
        self.outCoords = []
        self.outVals = []
        self.prevCoord = -1
        self.accumulatedValue = 0
        
        self.currentFiberSet = None
        self.partialOrFinal = None
        self.outputPartialOrFinal = None
        self.inputPartialOrFinal = None
        self.row = None
        self.JobID = None
        self.inputFibers = None
        self.Adata = None
        
        
        self.memoryWastedCycles = 0
        self.wastedCycles = 0
        
        
    def connect(self, scheduler: "Scheduler", fiberCache: "FiberCache"):
        self.scheduler = scheduler
        self.fiberCache = fiberCache
        self.nextFiberSet = self.scheduler.reassign()


        if not self.nextFiberSet:
            return
        # nextFiberSet[1] is a boolean representing reading B fibers or partial fibers, and nextFiberSet[4] is a list of fibers to read
        for x in range(len(self.nextFiberSet[1])):
            self.fiberCache.fetch(self.nextFiberSet[1][x], self.nextFiberSet[4][x])
        
        
        
    def switchAndReassign(self):  

        # save the output in the fiberCache using its jobID as a key if its a partial output, 
        # or save it in DRAM with its row value if its a final output.
        if self.prevCoord != -1:
            self.outCoords.append(self.prevCoord)
            self.outVals.append(self.accumulatedValue)
        
        self.prevCoord = -1
        self.accumulatedValue = 0
    
        self.fiberCache.save(self.outputPartialOrFinal,self.JobID, self.row, self.outCoords, self.outVals)
        self.outCoords = []
        self.outVals = []
        
        if self.currentFiberSet != None:
            # unlocks the fibers assigned to the now finished job
            self.fiberCache.unlockFibers(self.inputPartialOrFinal, self.inputFibers)
            
            # marks the finished job as completed
            self.scheduler.markCompleted(self.JobID)
            
            
        # switch Fiber sets
        self.currentFiberSet = self.nextFiberSet

        self.nextFiberSet = self.scheduler.reassign()

        if self.currentFiberSet != None:
            self.outputPartialOrFinal = self.currentFiberSet[0]
            self.inputPartialOrFinal = self.currentFiberSet[1]
            self.row = self.currentFiberSet[2]
            self.JobID = self.currentFiberSet[3]
            self.inputFibers = self.currentFiberSet[4]
            self.Adata = self.currentFiberSet[5]
            
            # checkIfFiberInCache checks if a fiber is in the cache. If it isn't, it makes sure that it is 
            # being fetched from DRAM. Normally, it should have been fetched, but its possible that it could be replaced 
            # in the short time it wasn't used
            t = False
            for i in range(len(self.inputFibers)):
                if self.fiberCache.checkIfFiberInCache(self.inputPartialOrFinal[i], self.inputFibers[i]):
                    # since the fibers are 100% guaranteed to be used basically instantly, we can lock them in so we don't replace them
                    # on accident.
                    self.fiberCache.lockFiber(self.inputPartialOrFinal[i], self.inputFibers[i])
                else:
                    t = True
            self.dataNotInCacheFlag = t

            if not t:       
                # now, we put all the values into our inputFibers queues
                for inp in range(len(self.inputFibers)):
                    if self.inputPartialOrFinal[inp]:
                        # setup B fibers data transfer from fiberCache into this input number
                        # this fills in both the coords[inp] and the values[inp]
                        self.fiberCache.read(self,inp,self.inputFibers[inp])
                    else:
                        # setup partial fibers data transfer from fiberCache into this input number
                        # this fills in both the coords[inp] and the values[inp]
                        self.fiberCache.consume(self,inp,self.inputFibers[inp])
        else:
            self.currentFiberSet = None
            self.partialOrFinal = None
            self.outputPartialOrFinal = None
            self.inputPartialOrFinal = None
            self.row = None
            self.JobID = None
            self.inputFibers = None
            self.Adata = None
            if self.nextFiberSet == None:
                self.endFlag = True
            return

        if self.nextFiberSet == None:
            return
        
        # fetch data for the next set of fibers
        for x in range(len(self.nextFiberSet[1])):
            self.fiberCache.fetch(self.nextFiberSet[1][x], self.nextFiberSet[4][x])    
        
   
    
    def cycle(self):
        
        """
        there isn't anything mentioned in the paper, other than how many banks FiberCache has, for its data transfer rate.
        I think the assumption is that as long as the fiber is in the cache, we don't worry about streaming it to the PEs
        
        1) If the output is running dry, stop compute
        2) have a set of queues to hold input data
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
            
        if self.endFlag:
            self.wastedCycles += 1
            return
        if self.dataNotInCacheFlag:
            t = False
            for i in range(len(self.inputFibers)):
                # checkIfFiberInCache checks if a fiber is in the cache. If it isn't, it makes sure that it is 
                if self.fiberCache.checkIfFiberInCache(self.inputPartialOrFinal[i], self.inputFibers[i]):
                    self.fiberCache.lockFiber(self.inputPartialOrFinal[i], self.inputFibers[i])
                else:
                    t = True
            if t:
                self.memoryWastedCycles += 1
                return
            # if the data is all in the fiberCache, we populate the input fifos with the data.
            # we don't worry about data transfer times or limited bandwidth between the fiberCache and the PEs
            self.dataNotInCacheFlag = t
            # If all fibers are in the cache at the same time, we lock them in place
            
            # now, we put all the values into our inputFibers queues
            for inp in range(len(self.inputFibers)):
                if self.inputPartialOrFinal[inp]:
                    # setup B fibers data transfer from fiberCache into this input number
                    self.fiberCache.read(self,inp,self.inputFibers[inp])
                else:
                    # setup partial fibers data transfer from fiberCache into this input number
                    self.fiberCache.consume(self,inp,self.inputFibers[inp])



        # check all the input fifos to see if they are all empty. If they are, we're done with this round of merging
        checkIfEmpty = True
        # find the minimum value in all the fifos (coordinate, way)
        minCoordWay = (None,None)
        
        
        for way,queue in enumerate(self.coordFifos):
            if queue:
                checkIfEmpty = False
                if minCoordWay[0] == None or minCoordWay[0] > queue[0]:
                    minCoordWay = (queue[0],way)
        
        if checkIfEmpty:
            self.switchAndReassign()
            return
        
        
        # remove the minimum coord from its queue
        self.coordFifos[minCoordWay[1]].popleft()
        
        # The value associated with this index
        value = self.valueFifos[minCoordWay[1]].popleft() * self.Adata[minCoordWay[1]]
        
        #print("SMALLEST COORD, WAY OF SMALLEST COORD, VAL OF SMALLEST COORD",minCoordWay, value)

        
        #print("COORDS", self.coordFifos)
        #print("VALUES", self.valueFifos)
        self.total += value
        
        if self.prevCoord != minCoordWay[0]:
            if self.prevCoord != -1:
                self.outCoords.append(self.prevCoord)
                self.outVals.append(self.accumulatedValue)
            self.prevCoord = minCoordWay[0]
            self.accumulatedValue = value
        else:
            self.accumulatedValue += value
    


        
    
