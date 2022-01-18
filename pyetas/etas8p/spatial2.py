# -*- coding: utf-8 -*-
# pyetas
# Copyright (C) 2021-2022 Salvatore Iacoletti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
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


