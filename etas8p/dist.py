# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""

import numpy as np
from numba import jit


@jit(nopython=True)
def dist(x1, y1, x2, y2):
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

@jit(nopython=True)
def dist2(x1, y1, x2, y2):
    return (x1 - x2)**2 + (y1 - y2)**2

def norm(x, dim):
    sumv = 0
    for i in range(dim):
        sumv += x[i]**2
    return np.sqrt(sumv)



#%%


if __name__ == "__main__":
    
    x1 = np.array([1., 1.])
    y1 = np.array([1., 2.])
    x2 = 3.
    y2 = 2.
    
    print(dist(x1, y1, x2, y2)) # test ok with c++
    print(dist2(x1, y1, x2, y2))
    print(norm(x1, x1.shape[0]))
    
    for i in range(0,10000):
        dist(x1, y1, x2, y2)
        
    for i in range(0,10000):
        dist2(x1, y1, x2, y2)
        
    for i in range(0,10000):
        norm(x1, x1.shape[0])
        