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
# import multiprocessing as mp

from pyetas.etas8p.dist import dist
from pyetas.etas8p.poly import polyinteg, pGauss, dGauss
import pyetas.etas8p.lambdaf0 as m0
import pyetas.etas8p.lambdaf1 as m1
import pyetas.etas8p.lambdaf2 as m2
import pyetas.etas8p.lambdaf3 as m3
import pyetas.etas8p.lambdaf4 as m4
import pyetas.etas8p.lambdaf5 as m5
import pyetas.etas8p.lambdaf6 as m6
import pyetas.etas8p.lambdaf6f as m6f
import pyetas.etas8p.lambdaf7 as m7
import pyetas.etas8p.lambdaf7f as m7f
from pyetas.etas8p.etasfit import calc_timedep_mc
# from pyrisk.utils.gardner_knopoff_window import GardnerKnopoffWindowOrig



#%% ***************************************************************************


def cdeclust(theta, rbwd, revents, rpoly, tperiod, mmin, model, pb_fix, voronoi):

    # extract events
    t_arr = np.array(revents['tt'])
    x_arr = np.array(revents['xx'])
    y_arr = np.array(revents['yy'])
    m_arr = np.array(revents['mm'])
    bk = np.array(revents['bkgd'])
    pb = np.array(revents['prob'])
    lam = np.array(revents['lambd'])
    flag = np.array(revents['flag'])
    faults = np.array(revents['fault'])
    fault_mode = not np.all(faults == None)
    # events = [revents[j] for j in revents.keys()]
    N = t_arr.shape[0]


    # extract polygon information
    px = rpoly['px']
    py = rpoly['py']
    npoly = len(py)

    # extract time period information
    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']

    # extract bandwidthes
    bwd = np.array(rbwd)

    # extract model paramters
    tht = [j[1] for j in theta.items()]

    # this loop was vectorized for speed
    for i in range(0, N):
        r0 = dist(x_arr[i], y_arr[i], x_arr, y_arr)
        dg = dGauss(np.array(r0), bwd)
        s = np.sum(pb * np.array(dg)) # np.sum(pb[flag==-2] * np.array(dg)[flag==-2])
        # s = 0
        # for j in range(0, N):
        #     # r0 = dist(x[i], y[i], x[j], y[j])
        #     # s += pb[j] * dGauss(r0, bwd[j])
        #     s += pb[j] * dg[j]
        bk[i] = s / (tlength - tstart2)
            
    revents['bkgd'] = bk

    if model == 7:
        mc = []
        for jj in range(0,N):
            mc.append(calc_timedep_mc(jj, t_arr, m_arr, mmin))
        mc = np.array(mc)

    # loop modified for speed (polyinteg and clambdaj not optimized)
    # temp = [None]*N
    # for i in range(0, N):
    #     temp[i] = res[i][0]
    #     lam[i] = res[i][1]
    temp = [None]*N
    for i in range(0, N):
        
        if voronoi is None:
            temp[i] = polyinteg(pGauss, np.array([bwd[i]]), npoly, px, py, x_arr[i], y_arr[i])
        else: # numerical integration with precomputed Voronoi tassellation
            points, areas = voronoi["points"], voronoi["areas"]
            r = dist(x_arr[i], y_arr[i], points[:,0], points[:,1])
            temp[i] = np.min([1., np.sum(pGauss(r, np.array([bwd[i]]))*areas)])
        
        if model == 0:
            lam[i] = m0.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
        elif model == 1:
            lam[i] = m1.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
        elif model == 2:
            lam[i] = m2.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
        elif model == 3:
            lam[i] = m3.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
        elif model == 4:
            lam[i] = m4.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
        elif model == 5:
            lam[i] = m5.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
        elif model == 6 and not fault_mode:
            # gk = GardnerKnopoffWindowOrig()
            # ta = gk.calc(m_arr+mmin)[1]*364.75
            ta = np.array([5*365]*m_arr.shape[0])
            lam[i] = m6.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk, ta)
        elif model == 6 and fault_mode: # import functions for model 6 with fault geometry
            ta = np.array([5*365]*m_arr.shape[0])
            lam[i] = m6f.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk, ta,
                                  faults, mmin)
        elif model == 7 and not fault_mode:
            ta = np.array([5*365]*m_arr.shape[0])
            lam[i] = m7.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk, ta,
                                 mc)            
        elif model == 7 and fault_mode:
            ta = np.array([5*365]*m_arr.shape[0])
            lam[i] = m7f.clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk, ta,
                                  faults, mmin, mc)
        else:
            raise Exception('problem with model module')
            
    s = np.sum(pb * np.array(temp)) # np.sum(pb[flag==-2] * np.array(temp)[flag==-2])
    if not pb_fix:
        revents['prob'] = list((tht[0]**2 * np.array(bk)) / np.array(lam))
    revents['lambd'] = lam
    
    # # # old loop
    # # s = 0
    # # for i in range(0, N):
    # #     w = [bwd[i]]
    # #     s += pb[i] * polyinteg(pGauss, w, npoly, px, py, x[i], y[i])
    # #     lam[i] = clambdaj(tht, i, t_arr, x_arr, y_arr, m_arr, bk)
    # #     revents['prob'][i] = (tht[0]**2 * bk[i]) / lam[i]
    # #     revents['lambd'][i] = lam[i]

    out = dict(revents=revents, integ0=s)
    return out


#%% ***************************************************************************


def decluster(theta, rbwd, revents, rpoly, tperiod, mmin, ndiv, model, pb_fix, voronoi):
    tht = {j[0]: np.sqrt(j[1]) for j in theta.items()}
    cbkg = cdeclust(tht, rbwd, revents, rpoly, tperiod, mmin, model, pb_fix, voronoi)
    return cbkg



#%%


# if __name__ == "__main__":
    
#     import time
#     from pyetas.etas8p.voronoi import get_voronoi
#     from myutils.utils_pickle import load_pickle, save_pickle
    
#     theta = load_pickle('test/param1')
#     rdata = load_pickle('test/rdata')
#     rbwd = load_pickle('test/rbwd')
#     revents = rdata['revents']
#     revents["fault"] = [None]*len(revents["mm"])
#     rpoly = rdata['rpoly']
#     tperiod = rdata['tperiod']
#     tht = {j[0]: math.sqrt(j[1]) for j in theta.items()}

#     ttt = time.time()
#     cbkg = cdeclust(tht, rbwd, revents, rpoly, tperiod, mmin=3.5,
#                     model=5, pb_fix=False, voronoi=None)
#     print(time.time()-ttt)
    
#     # R result 872.8452
#     print(cbkg['integ0'])
#     # R results #TODO
#     print(sum(cbkg['revents']['bkgd']), sum(cbkg['revents']['prob']), sum(cbkg['revents']['lambd']))
    
    
    
#     # integration done with Voronoi's diagrams
#     points, areas, _ = get_voronoi(rpoly, min_prec=0.005)
#     voronoi = {"points": points,
#                "areas": areas}    
#     ttt = time.time()
#     cbkg = cdeclust(tht, rbwd, revents, rpoly, tperiod, mmin=3.5,
#                     model=5, pb_fix=False, voronoi=voronoi)
#     print(time.time()-ttt)
    
#     # R result 872.8452
#     print(cbkg['integ0'])
#     print(sum(cbkg['revents']['bkgd']), sum(cbkg['revents']['prob']), sum(cbkg['revents']['lambd']))


    
#     r1 = 2.23606798 
#     w = np.array([0.05])
#     print(pGauss(r1, w)) # test ok with c++
    
    