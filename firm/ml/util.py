from __future__ import absolute_import, division, print_function

import sys
import os
from math import ceil
from datetime import datetime

import numpy as np


### Project Functions ###

def log(message):
    print('ML|PID-%s [%s]: %s' % (os.getpid(), datetime.now(), message))
    sys.stdout.flush()
    

### General ###

IDENTITY_FUNCTION = lambda x: x

    
### Collection Helper Functions ###
    
def is_iterable(object):
    try:
        iter(object)
        return True
    except TypeError:
        return False
        
def split_to_size(full_list, size):
    n_sublists = int(ceil(len(full_list) / size))
    return [full_list[(i * size):((i + 1) * size)] for i in range(n_sublists)]

    
### Numpy Helper Functions ###
    
def create_random_mask(size, n_true):
    mask = np.zeros(size).astype(bool)
    true_indices = np.random.choice(size, n_true, replace = False)
    mask[true_indices] = True
    return mask


### Multiprocessing Helper Functions ###

def distirbute_with_callback(thread_pool, func, args, callback):

    async_results = []

    for arg in args:
        async_result = thread_pool.apply_async(func, [arg], callback = callback)
        async_results.append(async_result)

    for async_result in async_results:
        
        async_result.wait()
        
        if not async_result.successful():
            log('Failed result: %s' % async_result.get())
