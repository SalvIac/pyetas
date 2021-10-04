# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""


import numpy as np

from pyrisk.etas.etas8p.dist import dist, dist2
from myutils.utils_pickle import load_pickle, save_pickle
from pyrisk.etas.etas8p.spatial5 import fr
from pyrisk.etas.etas8p.poly import polyinteg


# # spatial density function and its derivatives
# def fr(r, w):
#     gamma = w[0]
#     D = w[1]
#     q = w[2]
#     mag = w[3]
#     sig = D * np.exp(gamma * mag)
#     return (1 - np.power(1 + r**2 / sig, 1 - q)) / (2 * np.pi)


# double dgamma_fr(double r, double w[])
# {
#   double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
#   return (1 -q) * pow(1 + r * r / sig, -q) * mag * r * r / sig / (2 * M_PI);
# }

# double dD_fr(double r, double w[])
# {
#   double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
#   return (1 - q) * pow(1 + r * r / sig, -q) / D * r * r / sig / (2 * M_PI);
# }

# double dq_fr(double r, double w[])
# {
#   double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
#   return pow(1 + r * r / sig, 1 - q) * log(1 + r * r / sig) / (2 * M_PI);
# }



# vectorized version with numpy
# approximating the integral of a function on a polygon region
def polyintegXX_new(func, funcpara, px, py, cx, cy, ndiv):
    _sum = 0.
    for k in range(0, len(px) - 1):
    
        dxx = (px[k + 1] - px[k]) / ndiv
        dyy = (py[k + 1] - py[k]) / ndiv
        
        x1 = px[k] + dxx * np.array(range(0, ndiv))
        y1 = py[k] + dyy * np.array(range(0, ndiv))
        x2 = px[k] + dxx * np.array(range(1, ndiv+1))
        y2 = py[k] + dyy * np.array(range(1, ndiv+1))

        det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy)
        _id = np.ones(det.shape[0])
        _id[det < 0.] = -1.
        
        r1 = dist(x1, y1, cx, cy)
        r2 = dist(x2, y2, cx, cy)
        
        # dd = np.array( [dist2(x1[l], y1[l], x2[l], y2[l]) for l in range(0, ndiv)] )
        dd = dist2(x1[0], y1[0], x2[0], y2[0])
        
        theta = (r1**2 + r2**2 - dd)/(2 * r1 * r2)
        theta[np.abs(theta) > 1.] = 1. - 1.0e-10
        theta = np.arccos(theta)
        
        r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                  y1 + r1/(r1 + r2) * (y2 - y1), cx, cy)

        temp = _id * (func(r1, funcpara)/6 + func(r0, funcpara) * 2/3 +
                      func(r2, funcpara)/6) * theta
        
        _sum1 = np.sum(temp[ np.logical_and(np.abs(det) >= 1.0e-10, r1 + r2 > 1.0e-20) ])
        _sum += _sum1
    return _sum



# original code in c++
# approximating the integral of a function on a polygon region
def polyintegXX(func, funcpara, px, py, cx, cy, ndiv):
    _sum = 0
    for k in range(0, len(px) - 1):
    
        dxx = (px[k + 1] - px[k]) / ndiv
        dyy = (py[k + 1] - py[k]) / ndiv
        for l in range(0, ndiv):
            x1 = px[k] + dxx * l
            y1 = py[k] + dyy * l
            x2 = px[k] + dxx * (l + 1)
            y2 = py[k] + dyy * (l + 1)
            det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy)
            
            if abs(det) < 1.0e-10:
                continue
            
            _id = 1
            if (det < 0.):
                _id = -1
            
            r1 = dist(x1, y1, cx, cy)
            r2 = dist(x2, y2, cx, cy)
            theta = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2)
            if abs(theta) > 1:
                theta = 1 - 1.0e-10
            
            theta = np.arccos(theta)
            
            if r1 + r2 > 1.0e-20:            
                r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                          y1 + r1/(r1 + r2) * (y2 - y1), cx, cy)
                _sum += _id * (func(r1, funcpara)/6 + func(r0, funcpara) * 2/3 +
                        func(r2, funcpara)/6) * theta
    return _sum



#%% lambdatemporal


def lambdatemporal(t, fit, use_old=False):
    if not use_old:
        out = cxxlambdtemp(t, fit.param, fit.revents, fit.catalog.rpoly,
                           fit.catalog.rtperiod, fit.integ0, fit.ndiv)
    else:
        out = cxxlambdtemp_old(t, fit.param, fit.revents, fit.catalog.rpoly,
                               fit.catalog.rtperiod, fit.integ0, fit.ndiv)
    return np.array(out)



# vectorized version with numpy
# temporal intensity function: integrating over the spatial domain
def cxxlambdtemp(tg, theta, revents, rpoly, tperiod, integ0, ndiv):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    
    px = np.array(rpoly['px'])
    py = np.array(rpoly['py'])
    npoly = px.shape[0]
    
    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']    

    N = t.shape[0]
    ng = len(tg)
    sinteg = [None]*N
    out = [None]*ng
    
    temp = list()
    for i in range(0, N):
        w = np.array([gamma, D, q, m[i]])
        # temp.append(polyintegXX_new(fr, w, px, py, x[i], y[i], ndiv))
        temp.append(polyinteg(fr, w, npoly, px, py, x[i], y[i]))
    sinteg = A * np.exp(alpha * m) * np.array(temp)
    
    for j in range(0, ng):
        _sums = (p - 1.)/c * np.power(1. + (tg[j] - t[t < tg[j]])/c, - p) * sinteg[t < tg[j]]
        _sum = np.sum(_sums)
        out[j] = mu * integ0 /(tlength - tstart2) + _sum
    return out



# original code in c++
# temporal intensity function: integrating over the spatial domain
def cxxlambdtemp_old(tg, theta, revents, rpoly, tperiod, integ0, ndiv):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    
    px = rpoly['px']
    py = rpoly['py']

    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']    

    N = t.shape[0]
    ng = len(tg)
    sinteg = [None]*N
    out = [None]*ng
    
    for i in range(0, N):
        w = [gamma, D, q, m[i]]
        sinteg[i] = A * np.exp(alpha * m[i]) * \
                    polyintegXX(fr, w, px, py, x[i], y[i], ndiv)
    
    for j in range(0, ng):
        _sum = 0
        for i in range(0, N):
            if t[i] < tg[j]:
                _sum += (p - 1)/c * np.power(1 + (tg[j] - t[i])/c, - p) * sinteg[i]
        out[j] = mu * integ0 /(tlength - tstart2) + _sum
    return out
    
    

#%% timetransform


def timetransform(fit, use_old=False):
    if not use_old:
        out = cxxtimetrans(fit.param, fit.revents, fit.catalog.rpoly,
                          fit.catalog.rtperiod, fit.integ0, fit.ndiv)
    else:
        out = cxxtimetrans_old(fit.param, fit.revents, fit.catalog.rpoly,
                          fit.catalog.rtperiod, fit.integ0, fit.ndiv)
    return np.array(out)
  # cxxtimetrans(fit$param, obj$revents, obj$rpoly, obj$rtperiod,
  #             fit$integ0, fit$ndiv)



# vectorized version with numpy
def cxxtimetrans(theta, revents, rpoly, tperiod, integ0, ndiv):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    
    px = np.array(rpoly['px'])
    py = np.array(rpoly['py'])
    npoly = px.shape[0]
    
    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']

    N = t.shape[0]
    sinteg = [None]*N
    out = [None]*N
    
    temp = list()
    for i in range(0, N):
        w = np.array([gamma, D, q, m[i]])
        # temp.append(polyintegXX_new(fr, w, px, py, x[i], y[i], ndiv))
        temp.append(polyinteg(fr, w, npoly, px, py, x[i], y[i]))
    sinteg = A * np.exp(alpha * m) * np.array(temp)
    
    for j in range(0, N):
        _sums1 = (1 - np.power(1 + (t[j] - t[:j][t[:j] > tstart2])/c, 1 - p)) * \
                 sinteg[:j][t[:j] > tstart2]
        _sums2 = (np.power(1 + (tstart2 - t[:j][t[:j] <= tstart2])/c, 1 - p) -
                  np.power(1 + (t[j] - t[:j][t[:j] <= tstart2])/c, 1 - p)) * \
                  sinteg[:j][t[:j] <= tstart2]
        _sum = np.sum(_sums1) + np.sum(_sums2)
        out[j] = mu * integ0 * (t[j] - tstart2) / (tlength - tstart2) + _sum
    
    # for j in range(0, N):
    #     _sums1 = (1 - np.power(1 + (t[j] - t[:j])/c, 1 - p)) * sinteg[:j]
    #     _sums2 = (np.power(1 + (tstart2 - t[:j])/c, 1 - p) -
    #               np.power(1 + (t[j] - t[:j])/c, 1 - p)) * sinteg[:j]
    #     _sum = np.sum(_sums1[t[:j] > tstart2]) + np.sum(_sums2[t[:j] <= tstart2])
    #     out[j] = mu * integ0 * (t[j] - tstart2) / (tlength - tstart2) + _sum
    return out



# original code in c++
# transformed times: \tau_i
def cxxtimetrans_old(theta, revents, rpoly, tperiod, integ0, ndiv):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    
    px = rpoly['px']
    py = rpoly['py']

    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']

    N = t.shape[0]
    sinteg = [None]*N
    out = [None]*N
    
    for i in range(0, N):
        w = [gamma, D, q, m[i]]
        sinteg[i] = A * np.exp(alpha * m[i]) * \
                    polyintegXX(fr, w, px, py, x[i], y[i], ndiv)
    
    for j in range(0, N):
        _sum = 0
        for i in range(0, j):
            if t[i] > tstart2:
                _sum += (1 - np.power(1 + (t[j] - t[i])/c, 1 - p)) * sinteg[i]
            else:
                _sum += (np.power(1 + (tstart2 - t[i])/c, 1 - p) -
                         np.power(1 + (t[j] - t[i])/c, 1 - p)) * sinteg[i]
        out[j] = mu * integ0 * (t[j] - tstart2) / (tlength - tstart2) + _sum
    return out



#%% lambdaspatial


def lambdaspatial(x, y, fit, use_old=False):
    if not use_old:
        out = cxxlambspat(x, y, fit.param, fit.revents, fit.catalog.rpoly,
                          fit.catalog.rtperiod, fit.bwd)
    else:
        out = cxxlambspat_old(x, y, fit.param, fit.revents, fit.catalog.rpoly,
                          fit.catalog.rtperiod, fit.bwd)
    return np.array(out)
  # cxxlambspat(x, y, fit$param, obj$revents, obj$rpoly,
  #             obj$rtperiod, fit$bwd)



# vectorized version with numpt
# spatial intensity function: integrating over the temporal domain
def cxxlambspat(xg, yg, theta, revents, rpoly, tperiod, bwd):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    bk = np.array(revents['bkgd'])
    pb = np.array(revents['prob'])
    bwd = np.array(bwd)
    
    px = rpoly['px']
    py = rpoly['py']

    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']  

    N = t.shape[0]
    ng = len(xg)
    out = [None]*ng
    
    for j in range(0, ng):
        gint = np.zeros(N)
        gint[t > tstart2] = 1 - np.power(1 + (tlength - t[t > tstart2])/c, 1 - p)
        gint[t <= tstart2] = np.power(1 + (tstart2 - t[t <= tstart2])/c, 1 - p) - \
                             np.power(1 + (tlength - t[t <= tstart2])/c, 1 - p)
        r2 = dist2(x, y, xg[j], yg[j])
        sig = D * np.exp(gamma * m)
        _sum = np.sum(A * np.exp(alpha * m) * gint * (q - 1) / (sig * np.pi) * \
                      np.power(1 + r2 / sig, - q))
        s1 = np.sum(np.exp(-r2/(2 * bwd**2)) / (2 * np.pi * bwd**2))
        cum_s1 = np.cumsum(np.exp(-r2/(2 * bwd**2)) / (2 * np.pi * bwd**2))
        s2 = np.sum(pb * cum_s1)
        out[j] = _sum + mu * s2/(tlength - tstart2)
    return out



# original code in c++
# spatial intensity function: integrating over the temporal domain
def cxxlambspat_old(xg, yg, theta, revents, rpoly, tperiod, bwd):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    bk = np.array(revents['bkgd'])
    pb = np.array(revents['prob'])
    
    px = rpoly['px']
    py = rpoly['py']

    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']  

    N = t.shape[0]
    ng = len(xg)
    out = [None]*ng
    
    for j in range(0, ng):
        _sum = 0
        s1 = 0
        s2 = 0
        for i in range(0, N):
            if t[i] > tstart2:
                gint = 1 - np.power(1 + (tlength - t[i])/c, 1 - p)
            else:
                gint = np.power(1 + (tstart2 - t[i])/c, 1 - p) - \
                       np.power(1 + (tlength - t[i])/c, 1 - p)
            r2 = dist2(xg[j], yg[j], x[i], y[i])
            sig = D * np.exp(gamma * m[i])
            _sum += A * np.exp(alpha * m[i]) * gint * (q - 1) / (sig * np.pi) * \
                    np.power(1 + r2 / sig, - q)
            s1 += np.exp(-r2/(2 * bwd[i]**2)) / (2 * np.pi * bwd[i]**2)
            s2 += pb[i] *  s1
        out[j] = _sum + mu * s2/(tlength - tstart2)
    return out



#%%


if __name__ == "__main__":

    import time  
    
    etas_iran = load_pickle('test/etas_iran')
    fit = etas_iran


    #%% general data

    theta = fit.param
    revents = fit.revents
    rpoly = fit.catalog.rpoly
    tperiod = fit.catalog.rtperiod
    integ0 = fit.integ0
    ndiv = fit.ndiv
    bwd = fit.bwd
    
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    
    px = rpoly['px']
    py = rpoly['py']

    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = theta['mu']
    A = theta['A']
    c = theta['c']
    alpha = theta['alpha']
    p = theta['p']
    D = theta['D']
    q = theta['q']
    gamma = theta['gamma']

    N = t.shape[0]
    
    
    #%% test vectorized timetransform

    # t1 = time.time()
    # o11 = timetransform(fit, use_old=True)
    # print(time.time()-t1)
    
    # t2 = time.time()
    # o21 = timetransform(fit)
    # print(time.time()-t2)
    
    # print(np.allclose(o11, o21))
    

    #%% test vectorized lambdatemporal

    # tg = np.linspace(min(t), max(t), 1000)

    # t1 = time.time()
    # o12 = lambdatemporal(tg, fit, use_old=True)
    # print(time.time()-t1)

    # t2 = time.time()
    # o22 = lambdatemporal(tg, fit)
    # print(time.time()-t2)

    # print(np.allclose(o12, o22))


    #%% test vectorized lambdaspatial
    
    xg = np.linspace(-3.25, 30.7, 10)
    yg = np.linspace(-10.2, 5.7, 10)
    
    t1 = time.time()
    o1 = lambdaspatial(xg, yg, fit, use_old=True)
    print(time.time()-t1)

    t2 = time.time()
    o2 = lambdaspatial(xg, yg, fit)
    print(time.time()-t2)

    print(np.allclose(o1, o2))


    #%% test vectorized polyintegXX

    # for i in range(0, N):
    #     w = [gamma, D, q, m[i]]
    #     print(i, 100*(polyintegXX(fr, w, px, py, x[i], y[i], ndiv)/polyintegXX_new(fr, w, px, py, x[i], y[i], ndiv)-1))
    
    # ttt = time.time()
    # for i in range(0, N):
    #     w = [gamma, D, q, m[i]]
    #     polyintegXX_new(fr, w, px, py, x[i], y[i], ndiv)
    # print(time.time()-ttt)
    
    # ttt = time.time()
    # for i in range(0, N):
    #     w = [gamma, D, q, m[i]]
    #     polyintegXX(fr, w, px, py, x[i], y[i], ndiv)
    # print(time.time()-ttt)
    
    