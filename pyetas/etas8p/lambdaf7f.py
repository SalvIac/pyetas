# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""

import warnings

import numpy as np
from openquake.hazardlib.scalerel.wc1994 import WC1994

from pyrisk.etas.etas8p.dist import dist, dist2
from pyrisk.etas.etas8p.spatial5 import (fr, dq_fr, dgamma_fr, dD_fr,
                                         pdf_fr, dq_pdf_fr, dgamma_pdf_fr, dD_pdf_fr)
from pyrisk.etas.etas8p.poly import polyinteg
from pyrisk.etas.etas8p.lambdaf6 import (pdf_time_trunc, pdf_time_trunc_p, 
                                         pdf_time_trunc_c, integral_pdf_time_trunc,
                                         integral_pdf_time_trunc_p,
                                         integral_pdf_time_trunc_c)

                             
'''
this version of lambdaf6 is modified to include the fault geometry in the 
spatial PDF as in Guo et al. (2015) equation (9). it assumes that the fault is
discretized in a number of points and the productivity is evenly distributed
among these points. FYI: Guo et al. (2017) proposed a more advanced method but
it adds an additional variable. 
'''

#%% ***************************************************************************


def clambdaj(theta, j, t, x, y, m, bk, ta, faults=None, mmin=None, mc=None):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    # ta = 3*365 # 10yr
    
    part1 = A * np.exp(alpha * mc[:j]) # m[:j])
    
    # for loop modified for speed
    delta = t[j] - t[:j]
    sig = D * np.exp(gamma * m[:j])
    r2 = dist2(x[j], y[j], x[:j], y[:j])
    part_s = part1 * \
             pdf_time_trunc(delta, c, p, ta[:j]) * \
             (q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q)
    s = mu * bk[j] + np.sum(part_s)
    
    ################################ faults ###################################
    # the approach here is to modify s, deleting the point source contribution
    # and adding the contribution of all the points of the fault
    if faults is not None: # if faults vector provided
        for ii, fault in enumerate(faults[:j]): # loop through vector faults
            if fault is not None: # if a fault is found
            
                delta = t[j] - t[ii]
                part_kg = A * np.exp(alpha * mc[ii]) * \
                          pdf_time_trunc(delta, c, p, np.array([ta[ii]]))[0]
                          
                # delete contribution from previous loop
                sig = D * np.exp(gamma * m[ii])
                r2 = dist2(x[j], y[j], x[ii], y[ii])
                s = s -  part_kg * \
                         (q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q)
            
                # magnitude is scaled with Wells and Coppersmith (1994)
                num = fault['xx'].flatten().shape[0]
                wc = WC1994()
                mscaled = wc.get_median_mag(wc.get_median_area(m[ii]+mmin, None)/num, None) - mmin
                if mscaled < 0.:
                    warnings.warn("the scaled magnitude " +str(mscaled+mmin)+
                                  " of the fault points cannot be lower than Mmin: "+str(mmin))
                    mscaled = 0.
                
                # sum over all the points of the fault
                sig = D * np.exp(gamma * mscaled)
                r2 = dist2(x[j], y[j],
                           np.array(fault['xx'].flatten()), np.array(fault['yy'].flatten()))
                part_f = np.sum((q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q))
                
                # add contribution of each discretization point of fault
                s = s + 1/num * part_kg * part_f
    ################################ faults ###################################
    return s



#%% ***************************************************************************


def clambdajGr(theta, j, t, x, y, m, bk, fv, dfv, ta, faults=None, mmin=None,
               mc=None):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    # ta = 3*365 # 10yr
    
    # for loop modified for speed
    part1 = np.exp(alpha * mc[:j]) # m[:j])

    delta = t[j] - t[:j]
    part2 = pdf_time_trunc(delta, c, p, ta[:j])
    
    sig = D * np.exp(gamma * m[:j])
    r2 = dist2(x[j], y[j], x[:j], y[:j])
    part3 = (q - 1)/(sig * np.pi) * np.power(1 + r2/sig, - q)
    
    part_s = A * part1 * part2 * part3
    s = mu * bk[j] + np.sum(part_s)

    sg1 = bk[j]
    
    sg2 = np.sum(part1 * part2 * part3)
        
    part2_c = pdf_time_trunc_c(delta, c, p, ta[:j])
    sg3 = A * np.sum(part1 * part2_c * part3)
    
    part1_alpha = part1 * mc[:j]
    sg4 = A * np.sum(part1_alpha * part2 * part3)
    
    part2_p = pdf_time_trunc_p(delta, c, p, ta[:j])
    sg5 = A * np.sum(part1 * part2_p * part3)
    
    part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
    sg6 = A * np.sum(part1 * part2 * part3_d)
    
    part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
    sg7 = A * np.sum(part1 * part2 * part3_q)
    
    part3_gamma = part3 * (-m[:j] + q * m[:j] * (1 - 1/(1 + r2/sig)))
    sg8 = A * np.sum(part1 * part2 * part3_gamma)
    
    
    ################################ faults ###################################
    # the approach here is to modify s, deleting the point source contribution
    # and adding the contribution of all the points of the fault
    if faults is not None: # if faults vector provided
        for ii, fault in enumerate(faults[:j]): # loop through vector faults
            if fault is not None: # if a fault is found
            
                part1 = np.exp(alpha * mc[ii])
                
                delta = t[j] - t[ii]
                part2 = pdf_time_trunc(delta, c, p, np.array([ta[ii]]))[0]

                sig = D * np.exp(gamma * m[ii])
                r2 = dist2(x[j], y[j], x[ii], y[ii])
                part3 = (q - 1)/(sig * np.pi) * np.power(1 + r2/sig, - q)
                 
                part2_c = pdf_time_trunc_c(delta, c, p, np.array([ta[ii]]))[0]
                part1_alpha = part1 * mc[ii]
                part2_p = pdf_time_trunc_p(delta, c, p, np.array([ta[ii]]))[0]
                
                part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
                part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
                part3_gamma = part3 * (-m[ii] + q * m[ii] * (1 - 1/(1 + r2/sig)))

                # delete contribution from previous loop
                s = s - A * part1 * part2 * part3
                sg2 = sg2 - part1 * part2 * part3
                sg3 = sg3 - A * part1 * part2_c * part3
                sg4 = sg4 - A * part1_alpha * part2 * part3
                sg5 = sg5 - A * part1 * part2_p * part3
                sg6 = sg6 - A * part1 * part2 * part3_d
                sg7 = sg7 - A * part1 * part2 * part3_q
                sg8 = sg8 - A * part1 * part2 * part3_gamma

                # magnitude is scaled with Wells and Coppersmith (1994)
                num = fault['xx'].flatten().shape[0]
                wc = WC1994()
                mscaled = wc.get_median_mag(wc.get_median_area(m[ii]+mmin, None)/num, None) - mmin
                if mscaled < 0.:
                    warnings.warn("the scaled magnitude " +str(mscaled+mmin)+
                                  " of the fault points cannot be lower than Mmin: "+str(mmin))
                    mscaled = 0.
                
                # sum over all the points of the fault
                sig = D * np.exp(gamma * mscaled)
                r2 = dist2(x[j], y[j],
                           np.array(fault['xx'].flatten()), np.array(fault['yy'].flatten()))
                
                part3 = (q - 1) / (sig * np.pi) * np.power(1 + r2 / sig, - q)
                part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
                part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
                part3_gamma = part3 * (-mscaled + q * mscaled * (1 - 1/(1 + r2/sig)))                
                
                # add contribution of each discretization point of fault
                s = s + A * part1 * part2 * np.sum(1/num * part3)
                sg2 = sg2 + part1 * part2 * np.sum(1/num * part3)
                sg3 = sg3 + A * part1 * part2_c * np.sum(1/num * part3)
                sg4 = sg4 + A * part1_alpha * part2 * np.sum(1/num * part3)
                sg5 = sg5 + A * part1 * part2_p * np.sum(1/num * part3)
                sg6 = sg6 + A * part1 * part2 * np.sum(1/num * part3_d)
                sg7 = sg7 + A * part1 * part2 * np.sum(1/num * part3_q)
                sg8 = sg8 + A * part1 * part2 * np.sum(1/num * part3_gamma)
    ################################ faults ###################################
    
    fv = s
    dfv[ 0 ] = sg1 * 2 * theta[0]
    dfv[ 1 ] = sg2 * 2 * theta[1]
    dfv[ 2 ] = sg3 * 2 * theta[2]
    dfv[ 3 ] = sg4 * 2 * theta[3]
    dfv[ 4 ] = sg5 * 2 * theta[4]
    dfv[ 5 ] = sg6 * 2 * theta[5]
    dfv[ 6 ] = sg7 * 2 * theta[6]
    dfv[ 7 ] = sg8 * 2 * theta[7]    
    return fv, dfv



#%% ***************************************************************************


def cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, ta,
            faults=None, mmin=None, mc=None, voronoi=None):
    
    # if (m[j]+3.5)<=5.49:
    #     return 0.
    
    # extract model parameters
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    # ta = 3*365 # 10yr
    
    # time distribution integral
    sk = A * np.exp(alpha * mc[j]) #m[j])
    
    gi = integral_pdf_time_trunc(t[j], c, p, ta[j], tstart2, tlength)

    if faults is None: # no fault
        if voronoi is None:
            w = np.array([ gamma, D, q, m[j] ])
            si = polyinteg(fr, w, npoly, px, py, x[j], y[j])
        else: # numerical integration with precomputed Voronoi tassellation
            points, areas = voronoi["points"], voronoi["areas"]
            r = dist(x[j], y[j], points[:,0], points[:,1])
            si = np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])        
        
        
    

    if faults is not None: # if faults vector provided
        fault = faults[j]
        if fault is None: # if a fault is not found (same as before)
            if voronoi is None:
                w = np.array([ gamma, D, q, m[j] ])
                si = polyinteg(fr, w, npoly, px, py, x[j], y[j])
            else: # numerical integration with precomputed Voronoi tassellation
                points, areas = voronoi["points"], voronoi["areas"]
                r = dist(x[j], y[j], points[:,0], points[:,1])
                si = np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])  


    ################################ faults ###################################
        else: # if a fault is found
            # magnitude is scaled with Wells and Coppersmith (1994)
            num = fault['xx'].flatten().shape[0]
            wc = WC1994()
            mscaled = wc.get_median_mag(wc.get_median_area(m[j]+mmin, None)/num, None) - mmin
            if mscaled < 0.:
                warnings.warn("the scaled magnitude " +str(mscaled+mmin)+
                              " of the fault points cannot be lower than Mmin: "+str(mmin))
                mscaled = 0.
            
            # sum over all the points of the fault #TODO very slow
            si = 0.
            w = np.array([ gamma, D, q, mscaled ])
            for jj, (xx, yy, dd) in enumerate(zip(fault['xx'].flatten(),
                                                  fault['yy'].flatten(),
                                                  fault['depth'].flatten())):
                
                if voronoi is None:
                    si = si + polyinteg(fr, w, npoly, px, py, xx, yy)
                else: # numerical integration with precomputed Voronoi tassellation
                    raise Exception("not working")
                    # points, areas = voronoi["points"], voronoi["areas"]
                    # r = dist(xx, yy, points[:,0], points[:,1])
                    # si = si + np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])  

            si = 1/num * si
    ################################ faults ###################################

    return sk * gi * si



#%% ***************************************************************************


def cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv,
              ta, faults=None, mmin=None, mc=None, voronoi=None):

    # if (m[j]+3.5)<=5.49:
    #     return 0., [0.]*len(theta)

    # extract model parameters
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    # ta = 3*365 # 10yr
    
    sk = A * np.exp(alpha * mc[j]) #m[j])
    
    # time distribution integral and derivatives
    gi = integral_pdf_time_trunc(t[j], c, p, ta[j], tstart2, tlength)
    gip = integral_pdf_time_trunc_p(t[j], c, p, ta[j], tstart2, tlength)
    gic = integral_pdf_time_trunc_c(t[j], c, p, ta[j], tstart2, tlength)

    if faults is None: # no fault
        if voronoi is None:
            w = np.array([ gamma, D, q, m[j] ])
            si      = polyinteg(fr, w, npoly, px, py, x[j], y[j])
            sid     = polyinteg(dD_fr, w, npoly, px, py, x[j], y[j])
            siq     = polyinteg(dq_fr, w, npoly, px, py, x[j], y[j])
            sigamma = polyinteg(dgamma_fr, w, npoly, px, py, x[j], y[j])
        else: # numerical integration with precomputed Voronoi tassellation
            points, areas = voronoi["points"], voronoi["areas"]
            r = dist(x[j], y[j], points[:,0], points[:,1])
            si      = np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])
            sid     = np.sum(dD_pdf_fr(r, gamma, D, q, m[j])*areas)
            siq     = np.sum(dq_pdf_fr(r, gamma, D, q, m[j])*areas)
            sigamma = np.sum(dgamma_pdf_fr(r, gamma, D, q, m[j])*areas)




    if faults is not None: # if faults vector provided
        fault = faults[j]
        if fault is None: # if a fault is not found (same as before)
            if voronoi is None:
                w = np.array([ gamma, D, q, m[j] ])
                si      = polyinteg(fr, w, npoly, px, py, x[j], y[j])
                sid     = polyinteg(dD_fr, w, npoly, px, py, x[j], y[j])
                siq     = polyinteg(dq_fr, w, npoly, px, py, x[j], y[j])
                sigamma = polyinteg(dgamma_fr, w, npoly, px, py, x[j], y[j])
            else: # numerical integration with precomputed Voronoi tassellation
                points, areas = voronoi["points"], voronoi["areas"]
                r = dist(x[j], y[j], points[:,0], points[:,1])
                si      = np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])
                sid     = np.sum(dD_pdf_fr(r, gamma, D, q, m[j])*areas)
                siq     = np.sum(dq_pdf_fr(r, gamma, D, q, m[j])*areas)
                sigamma = np.sum(dgamma_pdf_fr(r, gamma, D, q, m[j])*areas)

    ################################ faults ###################################
        else: # if a fault is found
            # magnitude is scaled with Wells and Coppersmith (1994)
            num = fault['xx'].flatten().shape[0]
            wc = WC1994()
            mscaled = wc.get_median_mag(wc.get_median_area(m[j]+mmin, None)/num, None) - mmin
            if mscaled < 0.:
                warnings.warn("the scaled magnitude " +str(mscaled+mmin)+
                              " of the fault points cannot be lower than Mmin: "+str(mmin))
                mscaled = 0.
            
            # sum over all the points of the fault #TODO very slow
            si      = 0.
            sid     = 0.
            siq     = 0.
            sigamma = 0.
            w = np.array([ gamma, D, q, mscaled ])
            for jj, (xx, yy, dd) in enumerate(zip(fault['xx'].flatten(),
                                                  fault['yy'].flatten(),
                                                  fault['depth'].flatten())):     
                if voronoi is None:
                    si      = si      + polyinteg(fr, w, npoly, px, py, xx, yy)
                    sid     = sid     + polyinteg(dD_fr, w, npoly, px, py, xx, yy)
                    siq     = siq     + polyinteg(dq_fr, w, npoly, px, py, xx, yy)
                    sigamma = sigamma + polyinteg(dgamma_fr, w, npoly, px, py, xx, yy)
                else: # numerical integration with precomputed Voronoi tassellation
                    raise Exception("not working")
                    # points, areas = voronoi["points"], voronoi["areas"]
                    # r = dist(xx, yy, points[:,0], points[:,1])
                    # si      = si      + np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])
                    # sid     = sid     + np.sum(dD_pdf_fr(r, gamma, D, q, m[j])*areas)
                    # siq     = siq     + np.sum(dq_pdf_fr(r, gamma, D, q, m[j])*areas)
                    # sigamma = sigamma + np.sum(dgamma_pdf_fr(r, gamma, D, q, m[j])*areas)                
                
            si      = 1/num * si
            sid     = 1/num * sid            
            siq     = 1/num * siq
            sigamma = 1/num * sigamma            
    ################################ faults ################################### 
    
    fv = sk * gi * si
    dfv[0] = 0.
    dfv[1] = sk * gi  * si / A        * 2 * theta[1]
    dfv[2] = sk * gic * si            * 2 * theta[2]
    dfv[3] = sk * gi  * si * mc[j]    * 2 * theta[3] #mc[j] 
    dfv[4] = sk * gip * si            * 2 * theta[4]
    dfv[5] = sk * gi  * sid           * 2 * theta[5]
    dfv[6] = sk * gi  * siq           * 2 * theta[6]
    dfv[7] = sk * gi  * sigamma       * 2 * theta[7]
    return fv, dfv


