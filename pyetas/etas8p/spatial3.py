# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""



import numpy as np
from numba import jit


#%% ***************************************************************************


@jit(nopython=True)
def fr(r, w):
    D = w[0]
    q = w[1]
    return (1 - (1 + r**2 / D)**(1 - q)) / (2 * np.pi)

@jit(nopython=True)
def dD_fr(r, w):
    D = w[0]
    q = w[1]
    return (1 - q) * (1 + r**2 / D)**(-q) / D * r**2 / D / (2 * np.pi)

@jit(nopython=True)
def dq_fr(r, w):
    D = w[0]
    q = w[1]
    return (1 + r**2 / D)**(1 - q) * np.log(1 + r**2 / D) / (2 * np.pi)



