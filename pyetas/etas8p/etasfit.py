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

from pyetas.etas8p.dist import norm, dist
# import pyetas.etas8p.lambdaf0 as m0
import pyetas.etas8p.lambdaf1 as m1
import pyetas.etas8p.lambdaf2 as m2
import pyetas.etas8p.lambdaf3 as m3
import pyetas.etas8p.lambdaf4 as m4
import pyetas.etas8p.lambdaf5 as m5
import pyetas.etas8p.lambdaf6 as m6
import pyetas.etas8p.lambdaf6f as m6f
import pyetas.etas8p.lambdaf7 as m7
import pyetas.etas8p.lambdaf7f as m7f
import pyetas.etas8p.lambdaf6f_schoenberg as m6fs
# from pyrisk.utils.gardner_knopoff_window import GardnerKnopoffWindowOrig



    
# this would not be appropriate for entire catalogues (only for sequences)
# because it assumes that all the previous events are connected to the jth one
def calc_timedep_mc(j, t, m, mmin):
    mags = m[:j]+mmin
    delta = t[j] - t[:j]
    ind = (delta < 1.) & (mags >= 6.)
    if np.any(ind):
        # estimated incompleteness function for California (Helmstetter et al. 2006)
        # where time is the time (in days) after an earthquake with magnitude mag
        mc = (mags[ind])-4.5-0.75*np.log10(delta[ind])
        mc = np.clip(mc, mmin, m[j] + mmin)
        mc_max = np.max(mc)
        return m[j] + mmin - mc_max
    else:
        return m[j]
            
#     if complete_mode is not None:
#         # estimated incompleteness function for California (Helmstetter and Shaw, 2006)
#         # where time is the time (in days) after an earthquake with magnitude mag
#         time_main = complete_mode["time_main"]
#         mag_main = complete_mode["mag_main"]
#         time = self.date2day(dt, time_main, tz=tz) # time after mainshock
#         dmc = np.clip(mag_main-4.5-0.75*np.log10(time[time > 0.]),
#                       mag_threshold, mm+mag_threshold) - mag_threshold
#         mm.loc[time > 0.] += -dmc
#         if np.any(mm<0.):
#             raise Exception("error with complete_mode!")




#%% ***************************************************************************


def clinesearch(rdata, xOld, h, fv, verbose, ram, model_module, voronoi):

    const2 = 1.0e-16
    xNew = [None]*len(xOld)
    # ram1, ram2, ram3, fv1, fv2, fv3, a1, a2, a3, b1, b2

    if ram <= 1.0e-30:
        ram = 0.1

    hnorm = norm(h, len(xOld))
    if hnorm > 1:
        ram = ram/hnorm

    ram1 = 0.
    ram2 = ram
    fv1  = fv

    for i in range(0, len(xOld)):
        xNew[i] = xOld[i] + ram2 * h[i]
    fv2 = cloglkhd(xNew, rdata, verbose, model_module, voronoi)

    if fv2 > fv1:
        # goto stat50
        while fv2 > fv1:
            ram3 = ram2
            fv3 = fv2
            ram2 = ram3 * 0.1
            if ram2 * hnorm < const2:
                ram = 0.
                return ram
            for i in range(0, len(xOld)):
                xNew[i] = xOld[i] + ram2 * h[i]
            fv2 = cloglkhd(xNew, rdata, verbose, model_module, voronoi)
    else:
        # stat30:
        ram3 = ram2*2.0
        for i in range(0, len(xOld)):
            xNew[i] = xOld[i] + ram3 * h[i]
        fv3 = cloglkhd(xNew, rdata, verbose, model_module, voronoi)
        while not (fv3 > fv2):
            ram1 = ram2
            ram2 = ram3
            fv1 = fv2
            fv2 = fv3
            ram3 = ram2*2.0
            for i in range(0, len(xOld)):
                xNew[i] = xOld[i] + ram3 * h[i]
            fv3 = cloglkhd(xNew, rdata, verbose, model_module, voronoi)

    # stat70:
    a1 = (ram3 - ram2) * fv1
    a2 = (ram1 - ram3) * fv2
    a3 = (ram2 - ram1) * fv3
    b2 = (a1 + a2 + a3) * 2.
    b1 = a1 * (ram3 + ram2) + a2 * (ram1 + ram3) + a3 * (ram2 + ram1)
    if (b2 == 0.):
        ram = ram2
        return ram
    else:
        ram = b1 / b2
        for i in range(0, len(xOld)):
            xNew[i] = xOld[i] + ram*h[i]
        fv = cloglkhd(xNew, rdata, verbose, model_module, voronoi)
        if (ram > ram2):
            if (fv <= fv2):
                ram1 = ram2
                ram2 = ram
                fv1 = fv2
                fv2 = fv
                # goto stat200130
            else:
                ram3 = ram
                fv3 = fv
                # goto stat200130
        else:
            if (fv >= fv2):
                 ram1 = ram
                 fv1 = fv
                 # goto stat200130
            else:
                 ram3 = ram2
                 ram2 = ram
                 fv3 = fv2
                 fv2 = fv
                 # goto stat200130
        
    # stat200130:
    a1 = (ram3 - ram2)*fv1
    a2 = (ram1 - ram3)*fv2
    a3 = (ram2 - ram1)*fv3
    b2 = (a1 + a2 + a3)*2.0
    b1 = a1 * (ram3 + ram2) + a2 * (ram1 + ram3) + a3 * (ram2 + ram1)
    if (b2 == 0.):
        ram = ram2
        return ram
    else:
        ram = b1 /b2
        for i in range(0, len(xOld)):
            xNew[i] = xOld[i] + ram*h[i]
        fv = cloglkhd(xNew, rdata, verbose, model_module, voronoi)
        if (fv2 < fv):
            ram = ram2
        return ram


#%% ***************************************************************************
# log-likelihood function of the model

def cloglkhd(tht, rdata, verbose, model_module, voronoi):
    # extract events
    t = np.array(rdata['revents']['tt'])
    x = np.array(rdata['revents']['xx'])
    y = np.array(rdata['revents']['yy'])
    m = np.array(rdata['revents']['mm'])
    flag = rdata['revents']['flag']
    bk = rdata['revents']['bkgd']
    events = [rdata['revents'][j] for j in rdata['revents'].keys()]
    N = len(events[0])
    fault = np.array(rdata['revents']['fault'])
    mmin = rdata['mmin'] # for faults

    # extract polygon information
    px = rdata['rpoly']['px']
    py = rdata['rpoly']['py']
    npoly = len(py)

    # extract time period information
    tstart2 = rdata['tperiod']['study_start']
    tlength = rdata['tperiod']['study_end']

    # extract integral of spatial intensity over the obsevation window
    integ0 = rdata['integ0']

    if "lambdaf6" in model_module.__name__ or "lambdaf7" in model_module.__name__:
        # gk = GardnerKnopoffWindowOrig()
        # ta = gk.calc(m+mmin)[1]*364.75
        ta = np.array([5*365]*m.shape[0])

    if "lambdaf7" in model_module.__name__:
        # mc here is the corrected relative magnitude (not completeness mag)
        mc = list()
        for jj in range(0,N):
            mc.append(calc_timedep_mc(jj, t, m, mmin))
        mc = np.array(mc)

    fv1 = 0.
    fv2 = 0.

    for j in range(0, N):
        
        if flag[j] == 1:
            if "lambdaf6" in model_module.__name__:
                if "f" == model_module.__name__[-1]: # fault mode
                    s = model_module.clambdaj(tht, j, t, x, y, m, bk, ta,
                                              fault, mmin)
                else:
                    s = model_module.clambdaj(tht, j, t, x, y, m, bk, ta)

            elif "lambdaf7" in model_module.__name__:
                if "f" == model_module.__name__[-1]: # fault mode
                    s = model_module.clambdaj(tht, j, t, x, y, m, bk, ta,
                                              fault, mmin, mc)
                else:
                    s = model_module.clambdaj(tht, j, t, x, y, m, bk, ta, mc)

            else:
                s = model_module.clambdaj(tht, j, t, x, y, m, bk)
                
            if s > 1.0e-25:
                fv1 += np.log(s)
            else:
                fv1 -= 100.0
        
        if "lambdaf6" in model_module.__name__:
            if "f" == model_module.__name__[-1]: # fault mode
                fv2 += model_module.cintegj(tht, j, t, x, y, m, npoly, px, py, 
                                            tstart2, tlength, ta, fault, mmin,
                                            voronoi)
            else:
                fv2 += model_module.cintegj(tht, j, t, x, y, m, npoly, px, py,
                                            tstart2, tlength, ta, voronoi)

        elif "lambdaf7" in model_module.__name__:
            if "f" == model_module.__name__[-1]: # fault mode
                fv2 += model_module.cintegj(tht, j, t, x, y, m, npoly, px, py, 
                                            tstart2, tlength, ta, fault, mmin,
                                            mc, voronoi)
            else:
                fv2 += model_module.cintegj(tht, j, t, x, y, m, npoly, px, py,
                                            tstart2, tlength, ta, mc, voronoi)

        else:
            fv2 += model_module.cintegj(tht, j, t, x, y, m, npoly, px, py,
                                        tstart2, tlength)

    fv2 += tht[0]**2 * (integ0)
    fv = -fv1 + fv2

    if verbose == 1:
        print("Function Value = {:8.5f}  {:7.2f}  {:7.2f}".format(fv, -fv1, fv2))

    return fv


#%% ***************************************************************************


def cloglkhdGr(tht, rdata, verbose, fv, dfv, model_module, voronoi):
    
    # extract events
    t = np.array(rdata['revents']['tt'])
    x = np.array(rdata['revents']['xx'])
    y = np.array(rdata['revents']['yy'])
    m = np.array(rdata['revents']['mm'])
    flag = rdata['revents']['flag']
    bk = rdata['revents']['bkgd']
    events = [rdata['revents'][j] for j in rdata['revents'].keys()]
    N = len(events[0])
    fault = np.array(rdata['revents']['fault'])
    mmin = rdata['mmin'] # for faults
    
    # extract polygon information
    px = rdata['rpoly']['px']
    py = rdata['rpoly']['py']
    npoly = len(py)

    # extract time period information
    tstart2 = rdata['tperiod']['study_start']
    tlength = rdata['tperiod']['study_end']

    # extract integral of spatial intensity over the obsevation window
    integ0 = rdata['integ0']

    if "lambdaf6" in model_module.__name__ or "lambdaf7" in model_module.__name__:
        # gk = GardnerKnopoffWindowOrig()
        # ta = gk.calc(m+mmin)[1]*364.75
        ta = np.array([5*365]*m.shape[0])

    if "lambdaf7" in model_module.__name__:
        mc = list()
        for jj in range(0,N):
            mc.append(calc_timedep_mc(jj, t, m, mmin))
        mc = np.array(mc)


    fv1 = 0.
    fv2 = 0.
    df1 = [0.]*len(tht)
    df2 = [0.]*len(tht)
    fv1temp = None
    g1temp = [None]*len(tht)
    fv2temp = None
    g2temp = [None]*len(tht)

    for j in range(0, N):
                
        if flag[j] == 1:

            if "lambdaf6" in model_module.__name__:
                if "f" == model_module.__name__[-1]: # fault mode
                    fv1temp, g1temp = model_module.clambdajGr(tht, j, t, x, y,
                                                              m, bk, fv1temp,
                                                              g1temp, ta, fault,
                                                              mmin)
                else:
                    fv1temp, g1temp = model_module.clambdajGr(tht, j, t, x, y,
                                                              m, bk, fv1temp,
                                                              g1temp, ta)

            elif "lambdaf7" in model_module.__name__:
                if "f" == model_module.__name__[-1]: # fault mode
                    fv1temp, g1temp = model_module.clambdajGr(tht, j, t, x, y,
                                                              m, bk, fv1temp,
                                                              g1temp, ta, fault,
                                                              mmin, mc)
                else:
                    fv1temp, g1temp = model_module.clambdajGr(tht, j, t, x, y,
                                                              m, bk, fv1temp,
                                                              g1temp, ta, mc)

            else:
                fv1temp, g1temp = model_module.clambdajGr(tht, j, t, x, y, m,
                                                          bk, fv1temp, g1temp)         

            if fv1temp > 1.0e-25:
                fv1 += np.log(fv1temp)
            else:
                fv1 -= 100.0
                
            for i in range(0, len(tht)):
                if fv1temp != 0.:
                    df1[i] += g1temp[i] / fv1temp
                # else:
                #     print("Warning fv1temp=0")

        if "lambdaf6" in model_module.__name__:
            if "f" == model_module.__name__[-1]: # fault mode
                fv2temp, g2temp = model_module.cintegjGr(tht, j, t, x, y, m,
                                                         npoly, px, py,
                                                         tstart2, tlength,
                                                         fv2temp, g2temp, ta,
                                                         fault, mmin, voronoi)
            else:
                fv2temp, g2temp = model_module.cintegjGr(tht, j, t, x, y, m,
                                                         npoly, px, py,
                                                         tstart2, tlength,
                                                         fv2temp, g2temp,
                                                         ta, voronoi)
                
        elif "lambdaf7" in model_module.__name__:
            if "f" == model_module.__name__[-1]: # fault mode
                fv2temp, g2temp = model_module.cintegjGr(tht, j, t, x, y, m,
                                                         npoly, px, py,
                                                         tstart2, tlength,
                                                         fv2temp, g2temp, ta,
                                                         fault, mmin, mc,
                                                         voronoi)
            else:
                fv2temp, g2temp = model_module.cintegjGr(tht, j, t, x, y, m,
                                                         npoly, px, py,
                                                         tstart2, tlength,
                                                         fv2temp, g2temp,
                                                         ta, mc, voronoi)                
        else:
            fv2temp, g2temp = model_module.cintegjGr(tht, j, t, x, y, m,
                                                     npoly, px, py, 
                                                     tstart2, tlength,
                                                     fv2temp, g2temp)
        fv2 += fv2temp
        for i in range(0, len(tht)):
	        df2[i] += g2temp[i]
    
    fv2 += tht[0]**2 * integ0
    df2[0] = integ0 * tht[0] * 2

    fv = -fv1 + fv2
    
    for i in range(0, len(tht)):
        dfv[i] = -df1[i] + df2[i]
        
    if verbose == 1:
        print("Function Value = {:8.5f}  {:7.2f}  {:7.2f}".format(fv, -fv1, fv2))
        for i in range(0, len(tht)):
            print("Gradiant[{:d}] = {:8.5f}    theta[{:d}] = {:2.8f}".format(i + 1, dfv[i], i + 1, tht[i]))
    
    return fv, dfv


# *******************************************************************************


def cfit(theta, rdata, ihess, rverbose, model_module, voronoi):
    
    # extract model parameters
    tht = [j[1] for j in theta.items()]
    
    # extract verbose control
    verbose = rverbose
    
    if verbose == 1:
        print("start Davidon-Fletcher-Powell procedure ...")

    tau1 = 1.0e-6
    tau2 = 1.0e-6
    eps1 = 1.0e-6
    eps2 = 1.0e-6
    const1 = 1.0e-17

    ramda = 0.05
    
    # h = [[0.]*8]*8
    s = [0.]*len(tht)
    dx = [0.]*len(tht)
    g0 = [0.]*len(tht)
    g = [0.]*len(tht)
    dg = [None]*len(tht)
    wrk = [None]*len(tht)
    fv = None
    
    # Initial estimate of inverse of hessian matrix
    h = ihess.copy()
    
    fv, g = cloglkhdGr(tht, rdata, verbose, fv, g, model_module, voronoi)

    for it in range(1, 10):
        for ic in range(0, 8):
            
            if ic > 0 or it > 1:
                
                for i in range(0, len(tht)):
                    dg[i] = g[i] - g0[i]
                
                for i in range(0, len(tht)):
                    sumv = 0.0
                    for j in range(0, len(tht)):
                        sumv += dg[j] * h[i][j]
                    wrk[i] = sumv
                
                s1 = 0.0
                s2 = 0.0
                for i in range(0, len(tht)):
                    s1 += wrk[i] * dg[i]
                    s2 += dx[i] * dg[i]
                
                if s1 <= const1 or s2 <= const1:
                    fvP = -fv
                    aicP = 2 * (fv + len(tht))
                    if verbose == 1:
                        print("loglikelihood = {:8.5f}    AIC = {:8.5f}".format(-fv, 2 * (fv + len(tht))))
                    # for i in range(0, 8):
                    #     dfvP[i] = g[i]
                    #     estimP[i] = tht[i]
                    #     for j in range(0, 8):
                    #         hessP[i][j] = h[i][j]
                    dfvP = g.copy()
                    estimP = tht.copy()
                    hessP = h.copy()
                    for i in range(0, len(tht)):
                        if verbose == 1:
                            print("theta[{:d}] = {:2.8f}    gradient[{:d}] = {:8.4f}".format(
                                  i + 1, pow(tht[i], 2), i + 1, g[i]))
                    out = dict(estimate = estimP,
                                fvout = fvP,
                                dfvout = dfvP,
                                aic = aicP,
                                hess = hessP)
                    return out

                if (s1 <= s2):
                    # fletcher type correction
                    stem = s1 / s2 + 1.0
                    for i, h_prev in enumerate(h):
                        h[i] = [h_prev[j] - (dx[i]*wrk[j] + wrk[i]*dx[j] - dx[i]*dx[j]*stem) / s2 for j in range(0, len(tht))]
                    # for i in range(0, 8):
                    #     for j in range(i, 8):
                    #         h[i][j] -= (dx[i] * wrk[j] + wrk[i] * dx[j] - dx[i] * dx[j] * stem) / s2
                    #         h[j][i] = h[i][j]
                else:
                    # Update the inverse of hessian matrix
                    for i, h_prev in enumerate(h):
                        h[i] = [h_prev[j] + dx[i]*dx[j]/s2 - wrk[i]*wrk[j]/s1 for j in range(0, len(tht))]
                        # for j in range(i, 8):
                        #     # davidon-fletcher-powell type correction
                        #     h[i][j] = h[i][j] + dx[i]*dx[j]/s2 - wrk[i]*wrk[j]/s1
                        #     h[j][i] = h[j][i] + dx[i]*dx[j]/s2 - wrk[i]*wrk[j]/s1

            ss = 0.0
            for i in range(0, len(tht)):
                sumv = 0.0
                for j in range(0, len(tht)):
                    sumv += h[i][j] * g[j]
                ss += sumv**2
                s[i] = -sumv
            
            s1 = 0.0
            s2 = 0.0
            for i in range(0, len(tht)):
                s1 += s[i] * g[i]
                s2 += g[i]**2
                
            ds2 = np.sqrt(s2)
            gtem = abs(s1) / ds2
            if gtem <= tau1 and ds2 <= tau2:
                
                fvP = -fv
                aicP = 2 * (fv + len(tht))
                if verbose == 1:
                    print("loglikelihood = {:8.5f}\tAIC = {:8.5f}\n".format(-fv, 2 * (fv + len(tht))))
                # for i in range(0, 8):
                #     dfvP[i] = g[i]
                #     estimP[i] = tht[i]
                #     for j in range(0, 8):
                #         hessP[i][j] = h[i][j]
                dfvP = g.copy()
                estimP = tht.copy()
                hessP = h.copy()
                for i in range(0, len(tht)):
                    if verbose == 1:
                        print("theta[{:d}] = {:2.8f}\t gradient[{:d}] = {:8.4f}".format(
                              i + 1, tht[i]**2, i + 1, g[i]))
                out = dict(estimate = estimP,
                           fvout = fvP,
                           dfvout = dfvP,
                           aic = aicP,
                           hess = hessP)
                return out

            if s1 >= 0:
                for i in range(0, len(tht)):
                    h[i] = [1. if i==j else 0. for j in range(0, len(tht))]
                    # for j in range(0, 8):
                    #     h[i][j] = 0.0
                    # h[i][i] = 1.0
                    s[i] = -s[i]
            
            ed = fv
            # line search
            if verbose == 1:
                print("\nstart line search along the specified direction ...")
            ramda = clinesearch(rdata, tht, s, ed, verbose, ramda, model_module, voronoi)
            if verbose == 1:
                print("back to Davidon-Fletcher-Powell Procedure: zeta = {:f}".format(ramda))

            s1 = 0.0
            
            # al = tht[3] # save alpha
            # qq = tht[6] # save q
            
            for i in range(0, len(tht)):
                dx[i] = s[i] * ramda
                
                # dx[6] = 1.e-9 # when imposing q=1.5
                # dx[3] = 1.e-9 # when imposing alpha=beta
                
                s1 += dx[i] * dx[i]
                g0[i] = g[i]
                tht[i] += dx[i]
                
                # tht[6] = qq # when imposing q=1.5
                # tht[3] = al # when imposing alpha=beta
                
            fv0 = fv
            fv, g = cloglkhdGr(tht, rdata, verbose, fv, g, model_module, voronoi)

            s2 = 0.0
            for i in range(0, len(tht)):
                s2 += g[i]**2
            
            if np.sqrt(s2) > tau2:
                continue
            if fv0/fv - 1.0 < eps1 and np.sqrt(s1) < eps2:
                # estimP = [None]*8
                # dfvP = [None]*8
                # hessP = [[None]*8]*8
                fvP = -fv
                aicP = 2 * (fv + len(tht))
                if verbose == 1:
                    print("loglikelihood = {:8.5f}\tAIC = {:8.5f}\n".format(-fv, 2 * (fv + len(tht))))
                dfvP = g.copy()
                estimP = tht.copy()
                hessP = h.copy()
                for i in range(0, len(tht)):
                    # dfvP[i] = g[i]
                    # estimP[i] = tht[i]
                    # for j in range(0, len(tht)):
                    #     hessP[i][j] = h[i][j]
                    if verbose == 1:
                        print("theta[{:d}] = {:2.8f}\t gradient[{:d}] = {:8.4f}".format(
                              i + 1, tht[i]**2, i + 1, g[i]))
                out = dict(estimate = estimP,
                           fvout = fvP,
                           dfvout = dfvP,
                           aic = aicP,
                           hess = hessP)
                return out
    return 0


#%% ***************************************************************************


def etasfit(theta, revents, rpoly, tperiod, integ0, m0, ihess, verbose, ndiv,
            eps, model, modified_etas, voronoi):
    tht = {j[0]: np.sqrt(j[1]) for j in theta.items()}
    # if False: #cxxcode
    #     cfit_res = cxxfit(tht, revents, rpoly, tperiod, integ0, ihess,
    #                   int(ndiv), eps, bool(verbose), int(nthreads))
    fault_mode = not np.all([np.array(revents["fault"]) == None])
    rdata = dict(revents=revents, rpoly=rpoly, tperiod=tperiod, integ0=integ0, mmin=m0)
    if model == 0: # import functions for model 0
        print('imported model 0 - space-independent')
        model_module = m0
    elif model == 1: # functions for model 1
        model_module = m1
    elif model == 2: # functions for model 2
        model_module = m2
    elif model == 3: # functions for model 3
        model_module = m3
    elif model == 4: # functions for model 4
        model_module = m4
    elif model == 5: # functions for model 5
        model_module = m5
        
    ####### functions for model 6 #####
    elif model == 6 and not fault_mode and not modified_etas: 
        model_module = m6
    elif model == 6 and fault_mode and not modified_etas: # functions for model 6 with fault geometry
        model_module = m6f
    ####### functions for model 6 modified (lim magnitude) ##### #TODO not working now!
    elif model == 6 and not fault_mode and modified_etas:
        model_module = m6
    elif model == 6 and fault_mode and modified_etas: # functions for model 6 with fault geometry
        model_module = m6f
    
    # identical to model 6 with complete magnitudes (Helmstetter et al. 2006)
    elif model == 7 and not fault_mode:
        model_module = m7
    elif model == 7 and fault_mode:
        model_module = m7f

    # identical to model 6 with Schoenberg 2013 methodology
    elif model == 99 and fault_mode:
        model_module = m6fs

    else:
        raise Exception('invalid model')
    print(model_module.__name__)

    x0 = [j[1] for j in tht.items()]
    cfit_res = cfit(tht, rdata, ihess, int(verbose), model_module, voronoi)
    
    if cfit_res == 0:
        raise Exception("Maximum Likelihood optimization failed to converge.\nPlease try a better starting point.")

    H = np.array(cfit_res['hess'])
    tht = np.array(cfit_res['estimate'])
    avcov = (1/4 * np.matmul(np.matmul(np.diag(1/tht), H), np.diag(1/tht))).tolist()
    
    cfit_res['estimate'] = [t**2 for t in tht]
    cfit_res['avcov'] = avcov
    
    cfit_res['loglik'] = cfit_res.pop('fvout')
    cfit_res['gradient'] = cfit_res.pop('dfvout')
    cfit_res['ihessian'] = cfit_res.pop('hess')

    return cfit_res
         


#%%


# if __name__ == "__main__":

#     import time
#     from pyetas.etas8p.voronoi import get_voronoi
#     from myutils.utils_pickle import load_pickle, save_pickle
    
#     theta = load_pickle('test/param1')
#     rdata = load_pickle('test/rdata')
#     revents = rdata['revents']
#     rpoly = rdata['rpoly']
#     tperiod = rdata['tperiod']
#     integ0 = rdata['integ0']
#     ihess = load_pickle('test/ihess')
#     rverbose = verbose = 1

#     ndiv = 1000
#     eps=1e-06
#     # cxxcode=False
#     # nthreads=1  
#     m0 = 4.5
#     revents["fault"] = [None]*len(revents["mm"])
#     rdata['mmin'] = m0
#     tht = [np.sqrt(j[1]) for j in theta.items()]
    
    



    
    
#     # # voronoi tassellation
#     # time1 = time.time()
#     # points, areas, _ = get_voronoi(rpoly, min_prec=0.05)
#     # voronoi = {"points": points,
#     #             "areas": areas}
#     # print("create vonoroi", time.time()-time1)
    




#     # print('\ncloglkhd original')
#     # print(cloglkhd(tht, rdata, verbose, m6, voronoi=None))
#     # time1 = time.time()
#     # for i in range(10):
#     #     cloglkhd(tht, rdata, False, m6, voronoi=None)
#     # print(time.time()-time1)  


#     # print('\ncloglkhd Vonoroi')
#     # print(cloglkhd(tht, rdata, verbose, m6, voronoi))
#     # time1 = time.time()
#     # for i in range(10):
#     #     cloglkhd(tht, rdata, False, m6, voronoi)
#     # print(time.time()-time1)




#     # print('\ncloglkhdGr original')
#     # fv = None
#     # dfv = [0.]*8
#     # print(cloglkhdGr(tht, rdata, verbose, fv, dfv, m6, voronoi=None))
#     # time1 = time.time()
#     # for i in range(10):
#     #     cloglkhdGr(tht, rdata, False, fv, dfv, m6, voronoi=None)
#     # print(time.time()-time1)


#     # print('\ncloglkhdGr Voronoi')
#     # fv = None
#     # dfv = [0.]*8
#     # print(cloglkhdGr(tht, rdata, verbose, fv, dfv, m6, voronoi))
#     # time1 = time.time()
#     # for i in range(10):
#     #     cloglkhdGr(tht, rdata, False, fv, dfv, m6, voronoi)
#     # print(time.time()-time1)




#     # print('\nclinesearch original')
#     # fv = 54496.430296197024
#     # ramda = 0.05
#     # s = [0.]*8
#     # print(clinesearch(rdata, tht, s, fv, verbose, ramda, m6, voronoi=None))
#     # time1 = time.time()
#     # for i in range(10):
#     #     clinesearch(rdata, tht, s, fv, False, ramda, m6, voronoi=None)
#     # print(time.time()-time1)  



#     #########################################################################
#     # test corrected relative magnitude (with mag compl)

#     # t = np.array(rdata['revents']['tt'])
#     # m = np.array(rdata['revents']['mm'])
#     # events = [rdata['revents'][j] for j in rdata['revents'].keys()]
#     # N = len(events[0])
#     # mmin = rdata['mmin'] # for faults
#     # t[1] = t[0]+0.001 # otherwise it does not produce any difference
#     # m[0] = 7-mmin
    
#     # mc = list()
#     # for jj in range(0,N):
#     #     mc.append(calc_timedep_mc(jj, t, m, mmin))
#     # mc = np.array(mc)

#     # import matplotlib.pyplot as plt
#     # plt.figure()
#     # plt.plot(t,m)
#     # plt.plot(t,mc)
#     # plt.show()

    
    
#     #########################################################################
#     print('\ncinteg model 6 schoenberg et al. vs ogata 1998')

#     t = np.array(rdata['revents']['tt'])
#     x = np.array(rdata['revents']['xx'])
#     y = np.array(rdata['revents']['yy'])
#     m = np.array(rdata['revents']['mm'])
#     flag = rdata['revents']['flag']
#     bk = rdata['revents']['bkgd']
#     events = [rdata['revents'][j] for j in rdata['revents'].keys()]
#     N = len(events[0])
#     fault = np.array(rdata['revents']['fault'])
#     mmin = rdata['mmin'] # for faults

#     # extract polygon information
#     px = rdata['rpoly']['px']
#     py = rdata['rpoly']['py']
#     npoly = len(py)

#     # extract time period information
#     tstart2 = rdata['tperiod']['study_start']
#     tlength = rdata['tperiod']['study_end']

#     # extract integral of spatial intensity over the obsevation window
#     integ0 = rdata['integ0']

#     ta = np.array([5*365]*m.shape[0])
    
#     # classic ogata 1998    
#     fv2 = 0.
#     for j in range(0, N):
#         fv2 += m6f.cintegj(tht, j, t, x, y, m, npoly, px, py, 
#                            tstart2, tlength, ta, fault, mmin)

#     # schoenberg et al.
#     fv2s = 0.
#     for j in range(0, N):
#         fv2s += m6fs.cintegj(tht, j, t, x, y, m, npoly, px, py, 
#                              tstart2, tlength, ta, fault, mmin)


#     stop
    
#     #########################################################################
#     print('\netasfit minimize model 7')
#     time1 = time.time()
#     print(theta)
#     res = etasfit(theta, revents, rpoly, tperiod, integ0, m0, ihess, verbose,
#                   ndiv, eps, model=7, modified_etas=False, voronoi=None)
#     print(time.time()-time1)
#     print(res)





#     ##########################################################################
#     # print('\netasfit minimize model 6')
#     # time1 = time.time()
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0, m0, ihess, verbose,
#     #               ndiv, eps, model=6, modified_etas=False, voronoi=None) # voronoi=None
#     # print(time.time()-time1)
#     # print(res)
#     # # results for m0=4.5
#     # # results = [0.5550567658536777, 0.15916320085945337, 0.03191834194207985,
#     # #             2.733122131710429, 1.0892490763379579, 0.01613177945248924,
#     # #             2.3170471458051125, 0.03270311617539947]



#     ##########################################################################
#     print('\netasfit minimize model 6 with Schoenberg 2013')
#     from scipy import optimize
#     time1 = time.time()
#     print(theta)
#     x0 = [np.sqrt(j[1]) for j in theta.items()]
#     res = optimize.minimize(cloglkhd, x0, args=(rdata, verbose, m6fs, None),
#                             method='nelder-mead',
#                             options={'xatol': 1e-6, 'disp': True})
#     tht = [t**2 for t in res.x]
#     print(time.time()-time1)
#     print(tht)
#     # results for m0=4.5
#     # results = [0.5550567658536777, 0.15916320085945337, 0.03191834194207985,
#     #             2.733122131710429, 1.0892490763379579, 0.01613177945248924,
#     #             2.3170471458051125, 0.03270311617539947]



    

#     ##########################################################################
#     # print('\netasfit minimize model 5')
#     # time1 = time.time()
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0,
#     #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=5)
#     # print(time.time()-time1)
#     # print(res)
#     # results = [0.555266260648366, 0.1864583473888933, 0.047195657755165474,
#                 # 2.7050217320762324, 1.1547553549260672, 0.015949991076632867,
#                 # 2.317120770715824, 0.023512565102052223]




#     ##########################################################################
#     # print('\netasfit minimize model 4')
#     # time1 = time.time()
#     # del theta['gamma']
#     # ihess = ihess[:7,:7]
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0,
#     #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=4)
#     # print(res)
#     # print(time.time()-time1)




#     ##########################################################################
#     # print('\netasfit minimize model 3')
#     # time1 = time.time()
#     # del theta['gamma']
#     # ihess = ihess[:7,:7]
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0,
#     #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=3)
#     # print(res)
#     # print(time.time()-time1)




#     ##########################################################################
#     # print('\netasfit minimize model 2')
#     # time1 = time.time()
#     # del theta['gamma']
#     # del theta['q']
#     # ihess = ihess[:6,:6]
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0,
#     #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=2)
#     # print(res)
#     # print(time.time()-time1)




#     ##########################################################################
#     # print('\netasfit minimize model 1')
#     # time1 = time.time()
#     # del theta['gamma']
#     # del theta['q']
#     # ihess = ihess[:6,:6]
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0,
#     #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=1)
#     # print(res)
#     # print(time.time()-time1)




#     ##########################################################################
#     # print('\netasfit minimize model 0')
#     # time1 = time.time()
#     # del theta['gamma']
#     # del theta['q']
#     # del theta['D']
#     # ihess = ihess[:5,:5]
#     # print(theta)
#     # res = etasfit(theta, revents, rpoly, tperiod, integ0,
#     #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=0)
#     # print(res)
#     # print(time.time()-time1)
    
