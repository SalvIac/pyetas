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

from pyetas.etas8p.dist import dist2
from pyetas.etas8p.poly import polyinteg



#%% ***************************************************************************


def clambdaj(theta, j, t, x, y, m, bk):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    # for loop modified for speed
    delta = t[j] - t[:j]
    part_s = A * np.exp(alpha * m[:j]) * \
             (p - 1)/c * np.power(1. + delta/c, -p)
    s = mu * bk[j] + np.sum(part_s)
    return s



#%% ***************************************************************************


def clambdajGr(theta, j, t, x, y, m, bk, fv, dfv):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    
    # for loop modified for speed
    part1 = np.exp(alpha * m[:j])
    
    delta = t[j] - t[:j]
    part2 = (p - 1)/c * np.power(1 + delta / c, - p)
    
    part_s = A * part1 * part2
    s = mu * bk[j] + np.sum(part_s)

    sg1 = bk[j]
    
    sg2 = np.sum(part1 * part2)
        
    part2_c = part2 * (-1/c - p/(c + delta) + p/c)
    sg3 = A * np.sum(part1 * part2_c)
    
    part1_alpha = part1 * m[:j]
    sg4 = A * np.sum(part1_alpha * part2)
    
    part2_p = part2 * (1/(p - 1) - np.log(1 + delta/c))
    sg5 = A * np.sum(part1 * part2_p)
    
    fv = s
    dfv[ 0 ] = sg1 * 2 * theta[0]
    dfv[ 1 ] = sg2 * 2 * theta[1]
    dfv[ 2 ] = sg3 * 2 * theta[2]
    dfv[ 3 ] = sg4 * 2 * theta[3]
    dfv[ 4 ] = sg5 * 2 * theta[4]
    return fv, dfv



#%% ***************************************************************************


def cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength):

    # extract model parameters
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2

    if t[j] > tstart2:
        ttemp = tlength - t[j]
        gi  = 1 - (1 + ttemp/c)**(1 - p)
    else:
        ttemp1 = tstart2 - t[j]
        ttemp2 = tlength - t[j]
        
        gi1  = 1 - (1 + ttemp1/c)**(1 - p)
        gi2  = 1 - (1 + ttemp2/c)**(1 - p)
        gi  = gi2 - gi1

    sk = A * np.exp(alpha * m[j])
    return sk * gi



#%% ***************************************************************************


def cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv):

    # extract model parameters
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    
    if t[j] > tstart2:
        ttemp = tlength - t[j]
        
        gi  = 1 - (1 + ttemp/c)**(1 - p)
        gic = - (1 - gi) * (1 - p) * ( 1/(c + ttemp) - 1/c)
        gip = - (1 - gi) * (np.log(c) - np.log(c + ttemp))

    else:
        ttemp1 = tstart2 - t[j]
        ttemp2 = tlength - t[j]
        
        gi1  = 1 - (1 + ttemp1/c)**(1 - p)
        gi2  = 1 - (1 + ttemp2/c)**(1 - p)
        gic1 = - (1 - gi1) * (1 - p) * (1/(c + ttemp1) - 1/c)
        gic2 = - (1 - gi2) * (1 - p) * (1/(c + ttemp2) - 1/c)
        gip1 = - (1 - gi1) * (np.log(c) - np.log(c + ttemp1))
        gip2 = - (1 - gi2) * (np.log(c) - np.log(c + ttemp2))
        
        gi  = gi2 - gi1
        gic = gic2 - gic1
        gip = gip2 - gip1
    
    sk = A * np.exp(alpha * m[j])
    fv = sk * gi
    dfv[0] = 0
    dfv[1] = sk * gi / A        * 2 * theta[1]
    dfv[2] = sk * gic           * 2 * theta[2]
    dfv[3] = sk * gi * m[j]     * 2 * theta[3]
    dfv[4] = sk * gip           * 2 * theta[4]
    return fv, dfv

