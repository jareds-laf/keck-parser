# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 16:03:56 2022

@author: Jared
"""

import os
from multiprocessing import Process
from multiprocessing.pool import Pool

def task(arg=":-)"):
    print(f"Bello\n{arg}", flush=True)
    
def f(x):
    return x*x
    
if __name__ == '__main__':
    # Define a task to run in a new process
    # arr = "hi :)"
    # p = Process(target=task(arr))
    
    # Start the task in a new process
    # p.start()
    
    # Wait for task to complete
    # p.join()
    
    # Pools are used to automatically manage worker processes
        # Worker processes: processes that are kept around for reuse rather than creating and destroying a bunch of processes
    
    # Create a process pool
    # Like when opening files, you have to pool.close() at the end to release resources. You can use with statements, too! :D
    with Pool(processes=4, initializer=f, initargs=[2]) as pool:
    # processes: number of owkrers to create and manage within the pool
    # initializer: name of function you want to run
    # initargs: arguments to the specified function
    # maxtasksperchild: number of tasks each worker can execute
        # child processes often accumulate resources without releasing them, so this argument is a good failsafe
    
    
    
    # Issue tasks for execution
        results = pool.map(task, range(4))
        
        
    array = [1,2,3,4,5,6,7,8]
    step = 1
    test = [array[i:i + step] for i in range(0, len(array), step)]
    print(test)
