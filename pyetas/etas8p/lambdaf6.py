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

from pyetas.etas8p.dist import dist, dist2
from pyetas.etas8p.spatial5 import (fr, dq_fr, dgamma_fr, dD_fr,
                                         pdf_fr, dq_pdf_fr, dgamma_pdf_fr, dD_pdf_fr)
from pyetas.etas8p.poly import polyinteg



########################## pdf of the delta time #############################

# from etas.simulation_functions (not imported because it's not super elegant) (t is a numpy array)
def pdf_time_trunc(t, c, p, ta):
    if not np.isclose(p, 1., rtol=1e-06):
        pdf = (1-p)/(np.power(c+ta, 1-p) - np.power(c, 1-p)) * np.power(c+t, -p)
    else:
        pdf = np.power(np.log(c+ta)-np.log(c), -1) * np.power(c+t, -1)
    pdf[t > ta] = 0.
    return pdf


# derivative of trunc pdf time in p (t is a numpy array)
def pdf_time_trunc_p(t, c, p, ta):
    if not np.isclose(p, 1., rtol=1e-06):
        # https://www.wolframalpha.com/input/?i=derive+%281-p%29%2F%28%28c%2Ba%29%5E%281-p%29+-+c%5E%281-p%29%29+*+%28c%2Bt%29%5E%28-p%29+in+p
        den = (ta+c)**(1-p) - c**(1-p)
        const = c**(1-p)*np.log(c) - (ta+c)**(1-p)*np.log(ta+c)
        ddp = pdf_time_trunc(t, c, p, ta) * (- 1/(1-p) - const/den - np.log(c+t))
         #-v1/(v3-v2) - ((1-p)*v1*(v2*np.log(c) - v3*np.log(ta+c)))/(v3-v2)**2 - \
         #     ((1-p)*v1*np.log(c+t))/(v3-v2)
    else:
        # ddp = np.zeros(t.shape)
        # to estimate the derivative (it should be zero, but is seems incorrect)
        d1 = pdf_time_trunc_p(t, c, p-1e-5, ta)
        d2 = pdf_time_trunc_p(t, c, p+1e-5, ta)
        ddp = (d1+d2)/2
    return ddp


# derivative of trunc pdf time in c (t is a numpy array)
def pdf_time_trunc_c(t, c, p, ta):
    if not np.isclose(p, 1., rtol=1e-06):
        # https://www.wolframalpha.com/input/?i=derive+%281-p%29%2F%28%28c%2Ba%29%5E%281-p%29+-+c%5E%281-p%29%29+*+%28c%2Bt%29%5E%28-p%29+in+c
        den = (ta+c)**(1-p) - c**(1-p)
        const = (1-p)*((ta+c)**(-p) - c**(-p))
        ddc = pdf_time_trunc(t, c, p, ta) * (- p/(c+t) - const/den)
    else:
        # https://www.wolframalpha.com/input/?i=derive+%28log%28c%2Ba%29-log%28c%29%29%5E%28-1%29+*+%28c%2Bt%29%5E%28-1%29+in+c
        ddc = pdf_time_trunc(t, c, p, ta) * (ta/(c*(ta+c)*(np.log(ta+c) - np.log(c))) - 1)
    return ddc



############## integral and derivatives of times (not delta) #################


# integral from tstart to tend (t here is a single value (not a numpy array))
def integral_pdf_time_trunc(t, c, p, ta, tstart, tend):
    if t > tend:
        raise Exception("error in integral_pdf_time_trunc: t > tend")
    # derive delta of the inputs
    dtstart = tstart - t
    dtend = tend - t

    # adjustments for tmax (ta)
    if dtend >= ta and dtstart >= ta:
        return 0.
    if dtend > ta:
        dtend = ta

    if not np.isclose(p, 1., rtol=1e-06):
        # https://www.wolframalpha.com/input/?i=integrate+%281-p%29%2F%28%28c%2Ba%29%5E%281-p%29+-+c%5E%281-p%29%29+*+%28c%2Bt%29%5E%28-p%29+dt
        den = (ta+c)**(1-p) - c**(1-p)
        if t > tstart:
            gi = np.power(c+dtend, 1-p)/den - c**(1-p)/den
        else:
            gi1 = np.power(c+dtstart, 1-p)/den - c**(1-p)/den
            gi2 = np.power(c+dtend, 1-p)/den - c**(1-p)/den
            gi = gi2 - gi1
    else:
        # https://www.wolframalpha.com/input/?i=integrate+%28log%28c%2Ba%29+-+log%28c%29%29%5E%28-1%29+*+%28c%2Bt%29%5E%28-1%29+dt
        den = np.log(ta+c) - np.log(c)
        if t > tstart:
            gi = (np.log(c + dtend) - np.log(c)) / den
        else:
            gi1 = (np.log(c + dtstart) - np.log(c)) / den
            gi2 = (np.log(c + dtend) - np.log(c)) / den
            gi = gi2 - gi1
    return gi


# derivative in p of the integral from tstart to tend (t here is a single value (not a numpy array))
def integral_pdf_time_trunc_p(t, c, p, ta, tstart, tend):
    if t > tend:
        raise Exception("error in integral_pdf_time_trunc: t > tend")
    # derive delta of the inputs
    delta = 0.
    dtstart = tstart - t
    dtend = tend - t

    # # adjustments for tmax (ta)
    if dtend >= ta and dtstart >= ta:
        return 0.
    if dtend > ta:
        dtend = ta

    if not np.isclose(p, 1., rtol=1e-06):
        # # https://www.wolframalpha.com/input/?i=derivate+%28c+%2B+t%29%5E%281+-+p%29%2F%28-c%5E%281+-+p%29+%2B+%28a+%2B+c%29%5E%281+-+p%29%29+in+p
        # cst = 1 / ((ta + c)**(1-p) - c**(1-p))
        # cst2 = c**(1-p)*np.log(c) - (ta+c)**(1-p)*np.log(ta+c)
        # if t > tstart:
        #     gip = ( - (cst**2)*cst2*(np.power(c+dtend, 1-p)) - cst*(np.power(c+dtend, 1-p)*np.log(c+dtend)) ) - \
        #           ( - (cst**2)*cst2*(np.power(c+delta, 1-p)) - cst*(np.power(c+delta, 1-p)*np.log(c+delta)) )
        # else:
        #     gip1 = ( - (cst**2)*cst2*(np.power(c+dtstart, 1-p)) - cst*(np.power(c+dtstart, 1-p)*np.log(c+dtstart)) ) - \
        #            ( - (cst**2)*cst2*(np.power(c+delta, 1-p)) - cst*(np.power(c+delta, 1-p)*np.log(c+delta)) )
        #     gip2 = ( - (cst**2)*cst2*(np.power(c+dtend, 1-p)) - cst*(np.power(c+dtend, 1-p)*np.log(c+dtend)) ) - \
        #            ( - (cst**2)*cst2*(np.power(c+delta, 1-p)) - cst*(np.power(c+delta, 1-p)*np.log(c+delta)) )
        #     gip = gip2 - gip1
            
            
            
        
            
        # if t[j] > tstart2:
        #     ttemp = tlength - t[j]
            
        #     gi  = 1 - (1 + ttemp/c)**(1 - p)
        #     gic = - (1 - gi) * (1 - p) * ( 1/(c + ttemp) - 1/c)
        #     gip = - (1 - gi) * (np.log(c) - np.log(c + ttemp))
    
        # else:
        #     ttemp1 = tstart2 - t[j]
        #     ttemp2 = tlength - t[j]
            
        #     gi1  = 1 - (1 + ttemp1/c)**(1 - p)
        #     gi2  = 1 - (1 + ttemp2/c)**(1 - p)
        #     gic1 = - (1 - gi1) * (1 - p) * (1/(c + ttemp1) - 1/c)
        #     gic2 = - (1 - gi2) * (1 - p) * (1/(c + ttemp2) - 1/c)
        #     gip1 = - (1 - gi1) * (np.log(c) - np.log(c + ttemp1))
        #     gip2 = - (1 - gi2) * (np.log(c) - np.log(c + ttemp2))
            
        #     gi  = gi2 - gi1
        #     gic = gic2 - gic1
        #     gip = gip2 - gip1            
            
        # https://www.wolframalpha.com/input/?i=derivate+%28%28c+%2B+t%29%5E%281+-+p%29+-+%28c%29%5E%281+-+p%29%29%2F%28-c%5E%281+-+p%29+%2B+%28a+%2B+c%29%5E%281+-+p%29%29+in+p
        den = (ta + c)**(1-p) - c**(1-p)
        cst = c**(1-p)*np.log(c) - (ta+c)**(1-p)*np.log(ta+c)
        if t > tstart:
            gip = (c**(1-p)*np.log(c) - (c+dtend)**(1-p)*np.log(c+dtend))/den - (((c+dtend)**(1-p) - c**(1-p))*cst)/den**2
        else:
            gip1 = (c**(1-p)*np.log(c) - (c+dtstart)**(1-p)*np.log(c+dtstart))/den - (((c+dtstart)**(1-p) - c**(1-p))*cst)/den**2
            gip2 = (c**(1-p)*np.log(c) - (c+dtend)**(1-p)*np.log(c+dtend))/den - (((c+dtend)**(1-p) - c**(1-p))*cst)/den**2
            gip = gip2 - gip1
    else:
        # # it should be zero
        # gip = 0.
        # to estimate the derivative
        d1 = integral_pdf_time_trunc_p(t, c, p-1e-5, ta, tstart, tend)
        d2 = integral_pdf_time_trunc_p(t, c, p+1e-5, ta, tstart, tend)
        gip = (d1+d2)/2
    return gip



# derivative in c of the integral from tstart to tend (t here is a single value (not a numpy array))
def integral_pdf_time_trunc_c(t, c, p, ta, tstart, tend):
    if t > tend:
        raise Exception("error in integral_pdf_time_trunc: t > tend")
    # derive delta of the inputs
    delta = 0.
    dtstart = tstart - t
    dtend = tend - t

    # # adjustments for tmax (ta)
    if dtend >= ta and dtstart >= ta:
        return 0.
    if dtend > ta:
        dtend = ta

    if not np.isclose(p, 1., rtol=1e-06):
        # # https://www.wolframalpha.com/input/?i=derivate+%28c+%2B+t%29%5E%281+-+p%29%2F%28-c%5E%281+-+p%29+%2B+%28a+%2B+c%29%5E%281+-+p%29%29+in+c
        # cst = 1 / ((ta + c)**(1-p) - c**(1-p))
        # cst2 = (1-p)*(ta+c)**(-p) - (1-p)*c**(-p)
        # if t > tstart:
        #     gic = cst * ((1-p)*(c+dtend)**(-p)) - cst**2 * (cst2 * np.power(c+dtend, 1-p)) - \
        #           cst * ((1-p)*(c+delta)**(-p)) - cst**2 * (cst2 * np.power(c+delta, 1-p))
        # else:
        #     gic1 = cst * ((1-p)*(c+dtstart)**(-p)) - cst**2 * (cst2 * np.power(c+dtstart, 1-p)) - \
        #           cst * ((1-p)*(c+delta)**(-p)) - cst**2 * (cst2 * np.power(c+delta, 1-p))
        #     gic2 = cst * ((1-p)*(c+dtend)**(-p)) - cst**2 * (cst2 * np.power(c+dtend, 1-p)) - \
        #           cst * ((1-p)*(c+delta)**(-p)) - cst**2 * (cst2 * np.power(c+delta, 1-p))
        #     gic = gic2 - gic1
        
        
        # https://www.wolframalpha.com/input/?i=derivate+%28%28c+%2B+t%29%5E%281+-+p%29+-+%28c%29%5E%281+-+p%29%29%2F%28-c%5E%281+-+p%29+%2B+%28a+%2B+c%29%5E%281+-+p%29%29+in+c
        den = (ta + c)**(1-p) - c**(1-p)
        cst2 = (1-p)*(ta+c)**(-p) - (1-p)*c**(-p)
        if t > tstart:
            gic = ((1-p)*(c+dtend)**(-p) - (1-p)*c**(-p))/den - (cst2 * ((c+dtend)**(1-p) - c**(1-p)))/den**2
        else:
            gic1 = ((1-p)*(c+dtstart)**(-p) - (1-p)*c**(-p))/den - (cst2 * ((c+dtstart)**(1-p) - c**(1-p)))/den**2
            gic2 = ((1-p)*(c+dtend)**(-p) - (1-p)*c**(-p))/den - (cst2 * ((c+dtend)**(1-p) - c**(1-p)))/den**2
            gic = gic2 - gic1        
        
    else:
        # https://www.wolframalpha.com/input/?i=derivative+%28log%28c+%2B+t%29-log%28c%29%29%2F%28-log%28c%29+%2B+log%28a+%2B+c%29%29+in+c
        den = np.log(ta+c) - np.log(c)
        cst = 1/(ta+c) - 1/c
        if t > tstart:
            gic = (1/(c+dtend)-1/c)/den - cst*(np.log(c+dtend) - np.log(c))/den**2
        else:
            gic1 = (1/(c+dtstart)-1/c)/den - cst*(np.log(c+dtstart) - np.log(c))/den**2
            gic2 = (1/(c+dtend)-1/c)/den - cst*(np.log(c+dtend) - np.log(c))/den**2
            gic = gic2 - gic1
    return gic







##############################  checks  ###################################

def checks():
    import matplotlib.pyplot as plt
    ta = 1*365
    delta = np.arange(0.0005, 1.2*365, 0.001)
    c = 0.005
    p = 1.1
    t = 50
    tstart = 100
    tend = 200
    deltat = np.array([2])
    pp = np.arange(0.5, 1.5, 0.0001)
    cc = np.arange(0.0001, 0.05, 0.0001)

    # check function
    vals = pdf_time_trunc(delta, c, p, ta)
    plt.figure()
    plt.semilogy(delta, vals, label="function")
    plt.show()
    print("integral entire pdf:", np.sum(vals)*np.min(np.diff(delta)))

    # check derivative of function in p
    vals = list()
    valsdp = list()
    for ppp in pp:
        vals.append( pdf_time_trunc(deltat, c, ppp, ta) )
        valsdp.append( pdf_time_trunc_p(deltat, c, ppp, ta) )
    plt.figure()
    plt.plot(pp, vals, label="function")
    plt.plot(pp, valsdp, label="derivative")
    plt.legend()
    plt.show()

    # check derivative of function in c
    vals = list()
    valsdp = list()
    for ccc in cc:
        vals.append( pdf_time_trunc(deltat, ccc, p, ta) )
        valsdp.append( pdf_time_trunc_c(deltat, ccc, p, ta) )
    plt.figure()
    plt.semilogy(cc, vals, label="function")
    plt.semilogy(cc, valsdp, label="derivative")
    plt.legend()
    plt.show()
    
    # check integral functions
    vals = pdf_time_trunc(delta, c, p, ta)
    valsp1 = pdf_time_trunc(delta, c, 1, ta)
    # sanity checks
    print("integral over the entire ta period (this has to be equal to 1):", integral_pdf_time_trunc(t, c, p, ta, t, ta+t)) # 
    print("p != 1 and t < tstart integral between 100 and 200:", integral_pdf_time_trunc(t, c, p, ta, tstart, tend))
    print("p != 1 and t < tstart discrete integral curve between 100 and 200:", np.sum(vals[(delta >= tstart-t) & (delta <= tend-t)])*np.min(np.diff(delta)))
    print("p == 1 and t < tstart integral between 100 and 200:", integral_pdf_time_trunc(t, c, 1, ta, tstart, tend))
    print("p == 1 and t < tstart discrete integral curve between 100 and 200:", np.sum(valsp1[(delta >= tstart-t) & (delta <= tend-t)])*np.min(np.diff(delta)))
    print("p != 1 and t > tstart integral between 100 and 200:", integral_pdf_time_trunc((tstart+tend)/2, c, p, ta, tstart, tend))
    print("p != 1 and t > tstart discrete integral curve between 100 and 200:", np.sum(vals[(delta >= tstart-(tstart+tend)/2) & (delta <= tend-(tstart+tend)/2)])*np.min(np.diff(delta)))
    print("p == 1 and t > tstart integral between 100 and 200:", integral_pdf_time_trunc((tstart+tend)/2, c, 1, ta, tstart, tend))
    print("p == 1 and t > tstart discrete integral curve between 100 and 200:", np.sum(valsp1[(delta >= tstart-(tstart+tend)/2) & (delta <= tend-(tstart+tend)/2)])*np.min(np.diff(delta)))
    print("p != 1 and tstart > ta and tend > ta integral:", integral_pdf_time_trunc(t, c, p, ta, ta+100, ta+200))
    print("p != 1 and tend > ta integral:", integral_pdf_time_trunc(t, c, p, ta, tstart, ta+200))

    # plot of different cases that might happen with t tstart tend and ta
    plt.figure()
    plt.axvspan(tstart, tend, color='r', alpha=0.1)
    t0 = -300
    plt.axvspan(t0, t0, color='b', alpha=0.5)
    plt.plot(delta+t0, pdf_time_trunc(delta, c, p, ta), label="function tstart tend both > ta")
    print("tstart tend both > ta integral:", integral_pdf_time_trunc(t0, c, p, ta, tstart, tend))
    t01 = -200
    plt.axvspan(t01, t01, color='b', alpha=0.5)
    plt.plot(delta+t01, pdf_time_trunc(delta, c, p, ta), label="function only tend > ta")
    print("only tend > ta integral:", integral_pdf_time_trunc(t01, c, p, ta, tstart, tend))
    t1 = -50
    plt.axvspan(t1, t1, color='b', alpha=0.5)
    plt.plot(delta+t1, pdf_time_trunc(delta, c, p, ta), label="function normal case 1")
    print("normal case 1 integral:", integral_pdf_time_trunc(t1, c, p, ta, tstart, tend))
    t2 = 50
    plt.axvspan(t2, t2, color='b', alpha=0.5)
    plt.plot(delta+t2, pdf_time_trunc(delta, c, p, ta), label="function normal case 2")
    print("normal case 2 integral:", integral_pdf_time_trunc(t2, c, p, ta, tstart, tend))
    t3 = 150
    plt.axvspan(t3, t3, color='b', alpha=0.5)
    plt.plot(delta+t3, pdf_time_trunc(delta, c, p, ta), label="function t in the middle of tstart-tend")
    print("t in the middle of tstart-tend integral:", integral_pdf_time_trunc(t3, c, p, ta, tstart, tend))
    t4 = 250
    plt.axvspan(t4, t4, color='b', alpha=0.5)
    plt.plot(delta+t4, pdf_time_trunc(delta, c, p, ta), label="function t > tend")
    print("t > tend should give an error (not possible in ETAS calibration)")
    plt.ylim([0,0.01])
    plt.legend()
    plt.show()
    
    # figure derivative of the integral in p
    vals = list()
    valsdp = list()
    for ppp in pp:
        vals.append( integral_pdf_time_trunc(deltat[0], c, ppp, ta, tstart, tend) )
        valsdp.append( integral_pdf_time_trunc_p(deltat[0], c, ppp, ta, tstart, tend) )
    plt.figure()
    plt.plot(pp, vals, label="integral")
    plt.plot(pp, valsdp, label="derivative in p")
    plt.legend()
    plt.show()
    
    # figure derivative of the integral in c
    vals = list()
    valsdc = list()
    for ccc in cc:
        vals.append( pdf_time_trunc(deltat, ccc, p, ta) )
        valsdc.append( integral_pdf_time_trunc_c(deltat, ccc, p, ta, tstart, tend) )
    plt.figure()
    plt.semilogy(cc, vals, label="integral")
    plt.semilogy(cc, valsdc, label="derivative in c")
    plt.legend()
    plt.show()

    
    # try if functions work with ta being a numpy array
    taarr = np.array([3*365]*delta.shape[0])
    vals = pdf_time_trunc(delta, c, p, ta)
    
    
    
    
    

#%% ***************************************************************************


def clambdaj(theta, j, t, x, y, m, bk, ta):
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
    
    part1 = A * np.exp(alpha * m[:j])
    
    # part1[(m[:j]+3.5)<=5.49] = 0.
    
    # for loop modified for speed
    delta = t[j] - t[:j]
    sig = D * np.exp(gamma * m[:j])
    r2 = dist2(x[j], y[j], x[:j], y[:j])
    part_s = part1 * \
             pdf_time_trunc(delta, c, p, ta[:j]) * \
             (q - 1) / (sig * np.pi) * np.power(1. + r2 / sig, -q)
    s = mu * bk[j] + np.sum(part_s)
    return s



#%% ***************************************************************************


def clambdajGr(theta, j, t, x, y, m, bk, fv, dfv, ta):
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
    part1 = np.exp(alpha * m[:j])
    # part1[(m[:j]+3.5)<=5.49] = 0.


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
    
    part1_alpha = part1 * m[:j]
    sg4 = A * np.sum(part1_alpha * part2 * part3)
    
    part2_p = pdf_time_trunc_p(delta, c, p, ta[:j])
    sg5 = A * np.sum(part1 * part2_p * part3)
    
    part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
    sg6 = A * np.sum(part1 * part2 * part3_d)
    
    part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
    sg7 = A * np.sum(part1 * part2 * part3_q)
    
    part3_gamma = part3 * (-m[:j] + q * m[:j] * (1 - 1/(1 + r2/sig)))
    sg8 = A * np.sum(part1 * part2 * part3_gamma)
    
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
            voronoi=None):

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
    gi = integral_pdf_time_trunc(t[j], c, p, ta[j], tstart2, tlength)

    if voronoi is None:
        w = np.array([ gamma, D, q, m[j] ])
        si = polyinteg(fr, w, npoly, px, py, x[j], y[j])
    else: # numerical integration with precomputed Voronoi tassellation
        raise Exception("not working")
        # points, areas = voronoi[j]["points"], voronoi[j]["areas"]
        # r = dist(x[j], y[j], points[:,0], points[:,1])
        # si = np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])
        
    sk = A * np.exp(alpha * m[j])
    return sk * gi * si



#%% ***************************************************************************


def cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv,
              ta, voronoi=None):

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
    
    # time distribution integral and derivatives
    gi = integral_pdf_time_trunc(t[j], c, p, ta[j], tstart2, tlength)
    gip = integral_pdf_time_trunc_p(t[j], c, p, ta[j], tstart2, tlength)
    gic = integral_pdf_time_trunc_c(t[j], c, p, ta[j], tstart2, tlength)

    if voronoi is None:
        w = np.array([ gamma, D, q, m[j] ])
        si      = polyinteg(fr, w, npoly, px, py, x[j], y[j])
        sid     = polyinteg(dD_fr, w, npoly, px, py, x[j], y[j])
        siq     = polyinteg(dq_fr, w, npoly, px, py, x[j], y[j])
        sigamma = polyinteg(dgamma_fr, w, npoly, px, py, x[j], y[j])
    else: # numerical integration with precomputed Voronoi tassellation
        points, areas = voronoi[j]["points"], voronoi[j]["areas"]
        r = dist(x[j], y[j], points[:,0], points[:,1])
        si      = np.min([1., np.sum(pdf_fr(r, gamma, D, q, m[j])*areas)])
        sid     = np.sum(dD_pdf_fr(r, gamma, D, q, m[j])*areas)
        siq     = np.sum(dq_pdf_fr(r, gamma, D, q, m[j])*areas)
        sigamma = np.sum(dgamma_pdf_fr(r, gamma, D, q, m[j])*areas)
    
    sk = A * np.exp(alpha * m[j])
    fv = sk * gi * si
    dfv[0] = 0.
    dfv[1] = sk * gi  * si / A        * 2 * theta[1]
    dfv[2] = sk * gic * si            * 2 * theta[2]
    dfv[3] = sk * gi  * si * m[j]     * 2 * theta[3]
    dfv[4] = sk * gip * si            * 2 * theta[4]
    dfv[5] = sk * gi  * sid           * 2 * theta[5]
    dfv[6] = sk * gi  * siq           * 2 * theta[6]
    dfv[7] = sk * gi  * sigamma       * 2 * theta[7]
    return fv, dfv



#%% ***************************************************************************


if __name__ == "__main__":
    
    import time
    from pyetas.etas8p.voronoi import get_mixed_voronoi

    theta = [np.sqrt(j[1]) for j in {'mu': 0.11459513439367879,
                                    'A': 0.01,
                                    'c': 0.01,
                                    'alpha': 1,
                                    'p': 1.3,
                                    'D': 0.01,
                                    'q': 2,
                                    'gamma': 1}.items()]
    j = 6
    t = np.array([0.0,
                  2.748414351851852,
                  2.771030092592593,
                  6.741145833333333,
                  14.265266203703703,
                  17.982349537037038,
                  22.867060185185185,
                  26.31388888888889,
                  27.278854166666665,
                  27.652256944444446])
    x = np.array([4.883913500346288,
                  4.0649713749413,
                  4.292608153475305,
                  -2.1840749499622643,
                  -3.508624822356413,
                  -3.9438294419639286,
                  -3.822768432121458,
                  -5.685457151743186,
                  3.7841680967638913,
                  5.042522844792312])
    y = np.array([3.3641999999999967,
                  -0.13560000000000372,
                  0.3831999999999951,
                  -2.4313000000000002,
                  -2.642400000000002,
                  -4.836400000000001,
                  -2.4711,
                  -7.162600000000001,
                  -2.0771000000000015,
                  5.8470999999999975])
    m = np.array([0.09999999999999964,
                 1.0999999999999996,
                 0.7000000000000002,
                 0.0,
                 0.9000000000000004,
                 3.0,
                 3.20000000000000018,
                 0.20000000000000018,
                 0.7999999999999998,
                 1.9000000000000004])
    bk = np.array([0.05191824, 0.05189822, 0.05186122, 0.05190771, 0.05189338,
                    0.05188762, 0.05190069, 0.05190376, 0.05190912, 0.05189609])
    
    tstart2 = 0.0
    tlength = 27.652256944444446
    
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
    
    ta = 5*365*np.ones(m.shape)

    
    print('\nclambdaj')
    print(clambdaj(theta, j, t, x, y, m, bk, ta))

    print('\nclambdajGr')
    fv = None
    dfv = [None]*8
    print(clambdajGr(theta, j, t, x, y, m, bk, fv, dfv, ta))
    
    
    start_time = time.time()
    print('\ncintegj')
    print(cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, ta))
    
    print('\ncintegjGr')
    fv = None
    dfv = [None]*8
    print(cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv, ta))
    print(time.time()-start_time)



    # integration done with Voronoi's diagrams
    voronoi = list()
    for i, (xx, yy) in enumerate(zip(x,y)):
        points, areas, _ = get_mixed_voronoi({"px": px, "py": py}, [xx, yy], 
                                             min_prec=1, disc_num=100)
        voronoi.append({"points": points,
                        "areas": areas})
    
    start_time = time.time()
    print('\ncintegj Voronoi')
    print(cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, ta, voronoi))
    
    print('\ncintegjGr Voronoi')
    fv = None
    dfv = [None]*8
    print(cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv, ta, voronoi))    
    print(time.time()-start_time)
    
    