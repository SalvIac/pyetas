# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""


import numpy as np

from pyrisk.etas.etas8p.dist import dist2
from pyrisk.etas.etas8p.spatial3 import fr, dq_fr, dD_fr
from pyrisk.etas.etas8p.poly import polyinteg



#%% ***************************************************************************


def clambdaj(theta, j, t, x, y, m, bk, fault):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    # for loop modified for speed
    delta = t[j] - t[:j]
    sig = D
    if fault[j] is None:
        r2 = dist2(x[j], y[j], x[:j], y[:j])
        part_s = A * np.exp(alpha * m[:j]) * \
                  (p - 1)/c * np.power(1. + delta/c, -p) * \
                  (q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q)
    else:
        num = fault[j]['xx'].shape[0] * fault[j]['xx'].shape[1]
        part_s_list = list()
        for jj, (xx, yy, dd) in enumerate(zip(fault[j]['xx'].flatten(),
                                              fault[j]['yy'].flatten(),
                                              fault[j]['depth'].flatten())):
            r2 = dist2(xx, yy, x[:j], y[:j])
            part_s_list.append(A/num * np.exp(alpha * m[:j]) * \
                              (p - 1)/c * np.power(1. + delta/c, -p) * \
                              (q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q))

        part_s = np.concatenate(part_s_list)
    s = mu * bk[j] + np.sum(part_s)
    return s



#%% ***************************************************************************


def clambdajGr(theta, j, t, x, y, m, bk, fv, dfv, fault):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    
    # for loop modified for speed
    part1 = np.exp(alpha * m[:j])
    
    delta = t[j] - t[:j]
    part2 = (p - 1)/c * np.power(1 + delta / c, - p)

    sig = D
    
    sg1 = bk[j]
    dfv[ 0 ] = sg1 * 2 * theta[0]
        
    part2_c = part2 * (-1/c - p/(c + delta) + p/c)
    part1_alpha = part1 * m[:j]
    part2_p = part2 * (1/(p - 1) - np.log(1 + delta/c))
    
    if fault[j] is None:
        
        r2 = dist2(x[j], y[j], x[:j], y[:j])
        part3 = (q - 1)/(sig * np.pi) * np.power(1 + r2/sig, - q)
        
        part_s = A * part1 * part2 * part3
        s = mu * bk[j] + np.sum(part_s)
    
        sg2 = np.sum(part1 * part2 * part3)
        sg3 = A * np.sum(part1 * part2_c * part3)
        sg4 = A * np.sum(part1_alpha * part2 * part3)
        sg5 = A * np.sum(part1 * part2_p * part3)
        
        part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
        sg6 = A * np.sum(part1 * part2 * part3_d)
        
        part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
        sg7 = A * np.sum(part1 * part2 * part3_q)
    
        fv = s
        dfv[ 1 ] = sg2 * 2 * theta[1]
        dfv[ 2 ] = sg3 * 2 * theta[2]
        dfv[ 3 ] = sg4 * 2 * theta[3]
        dfv[ 4 ] = sg5 * 2 * theta[4]
        dfv[ 5 ] = sg6 * 2 * theta[5]
        dfv[ 6 ] = sg7 * 2 * theta[6]
    
    else:
        
        num = fault[j]['xx'].shape[0] * fault[j]['xx'].shape[1]
        s = list()
        sg2 = list() 
        sg3 = list()
        sg4 = list()
        sg5 = list()
        sg6 = list()
        sg7 = list()
        for jj, (xx, yy, dd) in enumerate(zip(fault[j]['xx'].flatten(),
                                              fault[j]['yy'].flatten(),
                                              fault[j]['depth'].flatten())):
    
            r2 = dist2(xx, yy, x[:j], y[:j])
            part3 = (q - 1)/(sig * np.pi) * np.power(1 + r2/sig, - q)
            
            part_s = A/num * part1 * part2 * part3
            s = np.sum(part_s)
            
            sg2.append( np.sum(part1 * part2 * part3) )
            sg3.append( A/num * np.sum(part1 * part2_c * part3) )
            sg4.append( A/num * np.sum(part1_alpha * part2 * part3) )
            sg5.append( A/num * np.sum(part1 * part2_p * part3) )
            
            part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
            sg6.append( A/num * np.sum(part1 * part2 * part3_d) )
            
            part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
            sg7.append( A/num * np.sum(part1 * part2 * part3_q) )
        
        fv = mu * bk[j] + np.sum(s)

        dfv[ 1 ] = np.sum(sg2) * 2 * theta[1]
        dfv[ 2 ] = np.sum(sg3) * 2 * theta[2]
        dfv[ 3 ] = np.sum(sg4) * 2 * theta[3]
        dfv[ 4 ] = np.sum(sg5) * 2 * theta[4]
        dfv[ 5 ] = np.sum(sg6) * 2 * theta[5]
        dfv[ 6 ] = np.sum(sg7) * 2 * theta[6]

    return fv, dfv





#%% ***************************************************************************


def cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fault):

    # extract model parameters
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    
    if t[j] > tstart2:
        ttemp = tlength - t[j]
        gi  = 1 - (1 + ttemp/c)**(1 - p)
    else:
        ttemp1 = tstart2 - t[j]
        ttemp2 = tlength - t[j]
        
        gi1  = 1 - (1 + ttemp1/c)**(1 - p)
        gi2  = 1 - (1 + ttemp2/c)**(1 - p)
        gi  = gi2 - gi1

    w = np.array([ D, q ])
    if fault[j] is None:
        sk = A * np.exp(alpha * m[j])
        si = polyinteg(fr, w, npoly, px, py, x[j], y[j])
        return sk * gi * si
    else:
        num = fault[j]['xx'].shape[0] * fault[j]['xx'].shape[1]
        out = list()
        for jj, (xx, yy, dd) in enumerate(zip(fault[j]['xx'].flatten(),
                                              fault[j]['yy'].flatten(),
                                              fault[j]['depth'].flatten())):
            sk = A/num * np.exp(alpha * m[j])
            si = polyinteg(fr, w, npoly, px, py, xx, yy)
            out.append(sk * gi * si)
        return np.sum(out)
        
        



#%% ***************************************************************************


def cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv, fault):

    # extract model parameters
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    
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

    w = np.array([ D, q ])
    if fault[j] is None:
        si      = polyinteg(fr, w, npoly, px, py, x[j], y[j])
        sid     = polyinteg(dD_fr, w, npoly, px, py, x[j], y[j])
        siq     = polyinteg(dq_fr, w, npoly, px, py, x[j], y[j])
        
        sk = A * np.exp(alpha * m[j])
        fv      = sk * gi * si
        dfv[0] = 0
        dfv[1] = sk * gi  * si / A        * 2 * theta[1]
        dfv[2] = sk * gic * si            * 2 * theta[2]
        dfv[3] = sk * gi  * si * m[j]     * 2 * theta[3]
        dfv[4] = sk * gip * si            * 2 * theta[4]
        dfv[5] = sk * gi  * sid           * 2 * theta[5]
        dfv[6] = sk * gi  * siq           * 2 * theta[6]
        return fv, dfv
    else:
        num = fault[j]['xx'].shape[0] * fault[j]['xx'].shape[1]
        fv = list()
        
        dfv1 = list()
        dfv2 = list()
        dfv3 = list()
        dfv4 = list()
        dfv5 = list()
        dfv6 = list()
        dfv7 = list()
        
        for jj, (xx, yy, dd) in enumerate(zip(fault[j]['xx'].flatten(),
                                              fault[j]['yy'].flatten(),
                                              fault[j]['depth'].flatten())):
            si      = polyinteg(fr, w, npoly, px, py, xx, yy)
            sid     = polyinteg(dD_fr, w, npoly, px, py, xx, yy)
            siq     = polyinteg(dq_fr, w, npoly, px, py, xx, yy)
            
            sk = A/num * np.exp(alpha * m[j])
            fv.append( sk * gi * si )
            dfv1.append(0)
            dfv2.append( sk * gi  * si / (A/num)  * 2 * theta[1] )
            dfv3.append( sk * gic * si            * 2 * theta[2] )
            dfv4.append( sk * gi  * si * m[j]     * 2 * theta[3] )
            dfv5.append( sk * gip * si            * 2 * theta[4] )
            dfv6.append( sk * gi  * sid           * 2 * theta[5] )
            dfv7.append( sk * gi  * siq           * 2 * theta[6] )
            
        return np.sum(fv), [np.sum(dfv1), np.sum(dfv2), np.sum(dfv3),
                            np.sum(dfv4), np.sum(dfv5), np.sum(dfv6),
                            np.sum(dfv7)]




