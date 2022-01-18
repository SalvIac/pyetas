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
from numba import jit, prange

from pyetas.etas8p.dist import dist

import pyetas.etas8p.spatial1 as sp1
import pyetas.etas8p.spatial2 as sp2
import pyetas.etas8p.spatial3 as sp3
import pyetas.etas8p.spatial5 as sp5



#%% ***************************************************************************


def frint_orig(func, funcpara, x1, y1, x2, y2, cx, cy):
    iid = 1
    det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy)
    if det < 0:
        iid = -1
    if abs(det) < 1.0e-10:
        return 0.
    r1 = dist(x1, y1, cx, cy)
    r2 = dist(x2, y2, cx, cy)
    r12 = dist(x1, y1, x2, y2)
    theta = (r1 * r1 + r2 * r2 - r12 * r12)/(2 * r1 * r2)
    if abs(theta) > 1:
        theta = 1. - 1.0e-10
    theta = np.arccos(theta)
    if r1 + r2 > 1.0e-20:
        x0 = x1 + r1/(r1 + r2) * (x2 - x1)
        y0 = y1 + r1/(r1 + r2) * (y2 - y1)
    else:
        return 0.
    r0 = dist(x0, y0, cx, cy)
    f1 = func(r1, funcpara)
    f2 = func(r0, funcpara)
    f3 = func(r2, funcpara)
    return iid * (f1/6 + (f2 * 2)/3 + f3/6) * theta


def frint(func, funcpara, x1, y1, x2, y2, cx, cy):
    out = np.repeat(np.nan, x1.shape[0])
    iid = np.repeat(1, x1.shape[0])
    det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy)
    
    iid[det < 0] = -1

    r1 = dist(x1, y1, cx, cy)
    r2 = dist(x2, y2, cx, cy)
    r12 = dist(x1, y1, x2, y2)
    theta = (r1 * r1 + r2 * r2 - r12 * r12)/(2 * r1 * r2)
    theta[abs(theta) > 1.] = 1 - 1.0e-10
    theta = np.arccos(theta)
    
    x0 = x1 + r1/(r1 + r2) * (x2 - x1)
    y0 = y1 + r1/(r1 + r2) * (y2 - y1)
    r0 = dist(x0, y0, cx, cy)
    f1 = func(r1, funcpara)
    f2 = func(r0, funcpara)
    f3 = func(r2, funcpara)

    out = iid * (f1/6 + (f2 * 2)/3 + f3/6) * theta
    out[r1 + r2 <= 1.0e-20] = 0.
    out[abs(det) < 1.0e-10] = 0.
    return out



#%% ***************************************************************************
# approximating the integral of a function on a polygon region


def polyinteg_orig(func, funcpara, npoly, px, py, cx, cy):
    ndiv = 1000
    sumv = 0.
    for j in range(0, npoly-1):
        dxx = (px[j + 1] - px[j]) / ndiv
        dyy = (py[j + 1] - py[j]) / ndiv
        for i in range(0, ndiv):
            x1 = px[j] + dxx * i
            y1 = py[j] + dyy * i
            x2 = px[j] + dxx * (i + 1)
            y2 = py[j] + dyy * (i + 1)
            sumv += frint_orig(func, funcpara, x1, y1, x2, y2, cx, cy)
    return sumv


def polyinteg_medium(func, funcpara, npoly, px, py, cx, cy): # medium fast
    ndiv = 1000
    sumv = 0.
    for j in range(0, npoly-1):
        dxx = (px[j + 1] - px[j]) / ndiv
        dyy = (py[j + 1] - py[j]) / ndiv
        x1 = px[j] + dxx * np.array(range(0, ndiv))
        y1 = py[j] + dyy * np.array(range(0, ndiv))
        x2 = px[j] + dxx * np.array(range(1, ndiv+1))
        y2 = py[j] + dyy * np.array(range(1, ndiv+1))
        sumv += np.sum(frint(func, funcpara, x1, y1, x2, y2, cx, cy))
    return sumv



#%% ***************************************************************************



# value of the gaussian kernel function at r with bandwidth sig
@jit(nopython=True)
def dGauss(r, sig):
    return np.exp(-(r**2) /(2 * sig**2)) / (2 * np.pi * sig**2)

# integral of the gaussian kernel with bandwidth w[0] from 0 to r
@jit(nopython=True)
def pGauss(r, w):
    return (1 - np.exp(-(r**2) / (2 * w[0]**2))) / (2 * np.pi)



#%% ***************************************************************************



def polyinteg(func, funcpara, npoly, px, py, cx, cy): # fast
    ndiv = 1000
    dxx = np.diff(px) / ndiv
    dyy = np.diff(py) / ndiv
    x1 = np.tile(px[0:-1], (ndiv,1)) + np.dot(np.array(range(0, ndiv)).reshape(ndiv, 1), dxx.reshape(1, dxx.shape[0]))
    y1 = np.tile(py[0:-1], (ndiv,1)) + np.dot(np.array(range(0, ndiv)).reshape(ndiv, 1), dyy.reshape(1, dyy.shape[0]))
    x2 = np.tile(px[0:-1], (ndiv,1)) + np.dot(np.array(range(1, ndiv+1)).reshape(ndiv, 1), dxx.reshape(1, dxx.shape[0]))
    y2 = np.tile(py[0:-1], (ndiv,1)) + np.dot(np.array(range(1, ndiv+1)).reshape(ndiv, 1), dyy.reshape(1, dyy.shape[0]))
    sumv = loop1(func.__name__, func.__module__, funcpara, npoly, x1, y1, x2, y2, cx, cy)
    # for j in range(0, npoly-1):
    #     sumv += np.sum(frint(func, funcpara, x1[:,j], y1[:,j], x2[:,j], y2[:,j], cx, cy))
    return sumv


@jit(nopython=True) # , parallel=True) # twice as fast with parallel
def loop1(funcname, funcmodule, funcpara, npoly, x1, y1, x2, y2, cx, cy):
    sumv = 0.
    for j in prange(0, npoly-1):
        sumv += np.sum(frint_numba(funcname, funcmodule, funcpara, x1[:,j], y1[:,j], x2[:,j], y2[:,j], cx, cy))
    return sumv


@jit(nopython=True)
def frint_numba(funcname, funcmodule, funcpara, x1, y1, x2, y2, cx, cy):
    out = np.repeat(np.nan, x1.shape[0])
    iid = np.repeat(1, x1.shape[0])
    det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy)
    
    iid[det < 0] = -1

    r1 = dist(x1, y1, cx, cy)
    r2 = dist(x2, y2, cx, cy)
    r12 = dist(x1, y1, x2, y2)
    theta = (r1 * r1 + r2 * r2 - r12 * r12)/(2 * r1 * r2)
    theta[np.abs(theta) > 1.] = 1 - 1.0e-10
    theta = np.arccos(theta)
    
    x0 = x1 + r1/(r1 + r2) * (x2 - x1)
    y0 = y1 + r1/(r1 + r2) * (y2 - y1)
    r0 = dist(x0, y0, cx, cy)

    if funcname == 'pGauss':
        f1 = pGauss(r1, funcpara)
        f2 = pGauss(r0, funcpara)
        f3 = pGauss(r2, funcpara)
    else:
        if 'spatial5' in funcmodule:
            if funcname == 'fr':
                f1 = sp5.fr(r1, funcpara)
                f2 = sp5.fr(r0, funcpara)
                f3 = sp5.fr(r2, funcpara)
            elif funcname == 'dq_fr':
                f1 = sp5.dq_fr(r1, funcpara)
                f2 = sp5.dq_fr(r0, funcpara)
                f3 = sp5.dq_fr(r2, funcpara)
            elif funcname == 'dgamma_fr':
                f1 = sp5.dgamma_fr(r1, funcpara)
                f2 = sp5.dgamma_fr(r0, funcpara)
                f3 = sp5.dgamma_fr(r2, funcpara)
            elif funcname == 'dD_fr':
                f1 = sp5.dD_fr(r1, funcpara)
                f2 = sp5.dD_fr(r0, funcpara)
                f3 = sp5.dD_fr(r2, funcpara)
        elif 'spatial3' in funcmodule:
            if funcname == 'fr':
                f1 = sp3.fr(r1, funcpara)
                f2 = sp3.fr(r0, funcpara)
                f3 = sp3.fr(r2, funcpara)
            elif funcname == 'dq_fr':
                f1 = sp3.dq_fr(r1, funcpara)
                f2 = sp3.dq_fr(r0, funcpara)
                f3 = sp3.dq_fr(r2, funcpara)
            elif funcname == 'dD_fr':
                f1 = sp3.dD_fr(r1, funcpara)
                f2 = sp3.dD_fr(r0, funcpara)
                f3 = sp3.dD_fr(r2, funcpara)
        elif 'spatial2' in funcmodule:
            if funcname == 'fr':
                f1 = sp2.fr(r1, funcpara)
                f2 = sp2.fr(r0, funcpara)
                f3 = sp2.fr(r2, funcpara)
            elif funcname == 'dgamma_fr':
                f1 = sp2.dgamma_fr(r1, funcpara)
                f2 = sp2.dgamma_fr(r0, funcpara)
                f3 = sp2.dgamma_fr(r2, funcpara)
            elif funcname == 'dD_fr':
                f1 = sp2.dD_fr(r1, funcpara)
                f2 = sp2.dD_fr(r0, funcpara)
                f3 = sp2.dD_fr(r2, funcpara)
        elif 'spatial1' in funcmodule:
            if funcname == 'fr':
                f1 = sp1.fr(r1, funcpara)
                f2 = sp1.fr(r0, funcpara)
                f3 = sp1.fr(r2, funcpara)
            elif funcname == 'dD_fr':
                f1 = sp1.dD_fr(r1, funcpara)
                f2 = sp1.dD_fr(r0, funcpara)
                f3 = sp1.dD_fr(r2, funcpara)
        else:
            raise Exception('problems with spatial component module')


    out = iid * (f1/6 + (f2 * 2)/3 + f3/6) * theta
    out[r1 + r2 <= 1.0e-20] = 0.
    out[np.abs(det) < 1.0e-10] = 0.
    return out



#%% ***************************************************************************


if __name__ == "__main__":
    
    px = [-7.0152516525115995,
          7.0152516525115765,
          7.0152516525115765,
          -7.0152516525115995,
          -7.0152516525115995]
    py = [-9.141648000000004,
          -9.141648000000004,
          9.141647999999996,
          9.141647999999996,
          -9.141648000000004]
    npoly = 5
    w = np.array([0.05])
    cx = 4.883913500346288
    cy = 3.3641999999999967
    
    # from pyetas.etas8p.decluster import pGauss
    # from pyetas.etas8p.lambdaf import fr, dq_fr, dD_fr, dgamma_fr 
    import time
    
    # print('\nfrint_orig')
    # x1 = -6.959129639291507
    # y1 = -9.141648000000004
    # x2 = -6.945099135986483
    # y2 = -9.141648000000004
    # print(frint_orig(pGauss, w, x1, y1, x2, y2, cx, cy)) # 9.418890233908686e-05 test ok with c++
    # time1 = time.time()
    # for i in range(100000):
    #     frint_orig(pGauss, w, x1, y1, x2, y2, cx, cy)
    # print(time.time()-time1)
    
    # print('\nfrint_new')
    # x1 = np.array([-6.959129639291507, -6.959129639291507])
    # y1 = np.array([-9.141648000000004, -9.141648000000004])
    # x2 = np.array([-6.945099135986483, -6.945099135986483])
    # y2 = np.array([-9.141648000000004, -9.141648000000004])
    # w = np.array([0.05, 1., 2., 3.])
    # print(frint(pGauss, w, x1, y1, x2, y2, cx, cy)) # 9.418890233908686e-05 test ok with c++
    # time1 = time.time()
    # for i in range(1000):
    #     frint(fr, w, x1, y1, x2, y2, cx, cy)
    # print(time.time()-time1)


    
    print('\npolyinteg_orig')
    print(polyinteg_orig(pGauss, w, npoly, px, py, cx, cy)) # 1.0000 test ok with c++
    time1 = time.time()
    for i in range(100):
        polyinteg_orig(pGauss, w, npoly, px, py, cx, cy)
    print(time.time()-time1)
    
    print('\npolyinteg_medium')
    print(polyinteg_medium(pGauss, w, npoly, px, py, cx, cy))
    time1 = time.time()
    for i in range(10000):
        polyinteg_medium(pGauss, w, npoly, px, py, cx, cy)
    print(time.time()-time1)  

    print('\npolyinteg_superfast')
    print(polyinteg(pGauss, w, npoly, np.array(px), np.array(py), cx, cy))
    time1 = time.time()
    for i in range(10000):
        polyinteg(pGauss, w, npoly, np.array(px), np.array(py), cx, cy)
    print(time.time()-time1)  
    
    
    # import math
    # import numpy as np
    # r = 0.1
    
    # time1 = time.time()
    # for i in range(100000):
    #     math.exp((r**2))
    # print(time.time()-time1)

    # time1 = time.time()
    # for i in range(100000):
    #     np.exp((r**2))
    # print(time.time()-time1)
    
    
