import numpy as np

class Scheduler:
    # Takes A rows from memory, and assigns them one by one to PEs
    def __init__(self) -> None:
        pass
    
    
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
    def __init__(self) -> None:
        pass
    
