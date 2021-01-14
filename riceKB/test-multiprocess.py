#!/usr/bin/env python

import pprint
from globalVars import *
from globalVars import base_vocab_ns
import re
import requests
import os, sys
import datetime
import json
import copy
import time
import requests
from multiprocessing import Process,Lock,Manager,active_children,Pool,cpu_count,Queue
# import multiprocessing as mp
import random

def rand_num(queue):
    num = random.random()
    queue.put(num)


def my_func(x):
    print mp.current_process()
    return x**x

def my_pool():
    pool = mp.Pool(mp.cpu_count()*2)
    result = pool.map(my_func,[4,2,3,5,3,2,1,2])
    result_set2 = pool.map(my_func, [4,6,5,4,6,3,23,4,6])

    print result
    print result_set2

if __name__ == "__main__":
    # main()
    queue = Queue()
    processes = [Process(target=rand_num,args=(queue,)) for x in range(4)]

    for p in processes:
        p.start()

    for p in processes:
        p.join()


    results = [queue.get() for p in processes]
    print results
# manager = Manager()
# self.result = manager.dict()
# if dbs == 'all':
#     name_db = ["rap", "msu7", "iric"]
# else:
#     name_db = dbs
#     self.result.setdefault('snpseek', manager.dict())
#      try:
#         p = Pool(processes=cpu_count()*2)
#         for i in range(len(dbs)):
#             p.apply_async(self.query,args=("snpseek","snpseek",[str(chro), str(start_pos), str(end_pos), name_db[i]],))
#     finally:
#         p.close()
#         p.join()