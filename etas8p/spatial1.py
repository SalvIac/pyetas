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
    return (1. - np.exp(-r**2/(2*D))) / (2 * np.pi)


@jit(nopython=True)
def dD_fr(r, w):
    D = w[0]
    return - r**2 * np.exp( -r**2/(2*D) )/ (4 * D**2 * np.pi)

