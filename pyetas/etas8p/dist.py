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
        