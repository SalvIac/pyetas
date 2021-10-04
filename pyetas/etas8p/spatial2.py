# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""



import numpy as np
from numba import jit


#%% ***************************************************************************


@jit(nopython=True)
def fr(r, w):
    alpha = w[0]
    D = w[1]
    mag = w[2]
    sig = D * np.exp(alpha * mag)
    return (1. - np.exp(-r**2/(2*sig))) / (2 * np.pi)


@jit(nopython=True)
def dD_fr(r, w):
    alpha = w[0]
    D = w[1]
    mag = w[2]
    sig = D * np.exp(alpha * mag)
    return - r**2 * np.exp(-r**2/(2*sig) - alpha * mag)/ (4 * D**2 * np.pi)


@jit(nopython=True)
def dgamma_fr(r, w): # I left dgamma as a name for simplicity
    alpha = w[0]
    D = w[1]
    mag = w[2]
    sig = D * np.exp(alpha * mag)
    return - mag * r**2 * np.exp(-r**2/(2*sig) - alpha*mag) / (4*np.pi*D)


