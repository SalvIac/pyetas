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

import time

import os
import numpy as np



def pGauss(r, w):
    return (1 - np.exp(-(r**2) / (2 * w[0]**2))) / (2 * np.pi)

def dist(x1, y1, x2, y2):
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

def dist2(x1, y1, x2, y2):
    return (x1 - x2)**2 + (y1 - y2)**2

def frint(func, funcpara, x1, y1, x2, y2, cx, cy):
    iid = 1
    det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy)
    if det < 0:
        iid = -1
    if abs(det) < 1.0e-10:
        return 0
    r1 = dist(x1, y1, cx, cy)
    r2 = dist(x2, y2, cx, cy)
    r12 = dist(x1, y1, x2, y2)
    theta = (r1 * r1 + r2 * r2 - r12 * r12)/(2 * r1 * r2)
    if abs(theta) > 1:
        theta = 1 - 1.0e-10
    theta = np.arccos(theta)
    if r1 + r2 > 1.0e-20:
        x0 = x1 + r1/(r1 + r2) * (x2 - x1)
        y0 = y1 + r1/(r1 + r2) * (y2 - y1)
    else:
        return 0
    r0 = dist(x0, y0, cx, cy)
    f1 = func(r1, funcpara)
    f2 = func(r0, funcpara)
    f3 = func(r2, funcpara)
    return iid * (f1/6 + (f2 * 2)/3 + f3/6) * theta


def clambdaj(theta, j, t, x, y, m, bk):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    # for loop modified for speed
    delta = t[j] - t[:j]
    sig = D * np.exp(gamma * m[:j])
    r2 = dist2(x[j], y[j], x[:j], y[:j])
    part_s = A * np.exp(alpha * m[:j]) * \
              (p - 1)/c * np.power(1. + delta/c, -p) * \
              (q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q)
    s = mu * bk[j] + np.sum(part_s)
    # # old loop
    # s = mu * bk[j]
    # vect_r2 = dist2(x[j], y[j], x, y)
    # for i in range(0, j):
    #     part1 = np.exp(alpha * m[i])
    #     delta = t[j] - t[i]
    #     part2 = (p - 1)/c * (1 + delta/c)**(-p)
    #     sig = D * np.exp(gamma * m[i])
    #     r2 = vect_r2[i]
    #     part3 = (q - 1) / (sig * np.pi) * (1 + r2 / sig)**(-q)
    #     s += A * part1 * part2 * part3
    return s



def polyinteg(func, funcpara, npoly, px, py, cx, cy):
    ndiv = 1000
    sumv = 0
    for j in range(0, npoly-1):
        dxx = (px[j + 1] - px[j]) / ndiv
        dyy = (py[j + 1] - py[j]) / ndiv
        for i in range(0, ndiv):
            x1 = px[j] + dxx * i
            y1 = py[j] + dyy * i
            x2 = px[j] + dxx * (i + 1)
            y2 = py[j] + dyy * (i + 1)
        sumv += frint(func, funcpara, x1, y1, x2, y2, cx, cy)
    return sumv



#%%

# if __name__=="__main__":
    
#     from myutils.utils_pickle import load_pickle, save_pickle
    
#     revents = load_pickle('events.pkl')
#     rpoly = load_pickle('rpoly.pkl')
#     rbwd = load_pickle('rbwd.pkl')

#     t = revents['tt']
#     x = revents['xx']
#     y = revents['yy']
#     m = revents['mm']
#     bk = revents['bkgd']
#     pb = revents['prob']
#     lam = revents['lambd']
#     N = len(t)
    
#     bwd = np.array(rbwd)
    
#     theta = {'mu': 0.11459513439367879, 'A': 0.01, 'c': 0.01, 'alpha': 1, 'p': 1.3, 'D': 0.01, 'q': 2, 'gamma': 1}
#     tht = [np.sqrt(j[1]) for j in theta.items()]
    

#     # extract polygon information
#     px = rpoly['px']
#     py = rpoly['py']
#     npoly = len(py)

#     s = 0
#     t_arr = np.array(t)
#     x_arr = np.array(x)
#     y_arr = np.array(y)
#     m_arr = np.array(m)
    
#     for i in range(0, 1000):
#         w = [bwd[i]]
#         s += pb[i] * polyinteg(pGauss, w, npoly, px, py, x[i], y[i])
#         lam[i] = clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
    
#     print(s)
#     print(sum(lam))   



