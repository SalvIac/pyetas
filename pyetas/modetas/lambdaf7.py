# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""


import numpy as np

from pyetas.etas8p.dist import dist, dist2
from pyetas.etas8p.spatial5 import (fr, dq_fr, dgamma_fr, dD_fr,
                                         pdf_fr, dq_pdf_fr, dgamma_pdf_fr, dD_pdf_fr)
from pyetas.etas8p.poly import polyinteg
from pyetas.etas8p.lambdaf6 import (pdf_time_trunc, pdf_time_trunc_p, 
                                         pdf_time_trunc_c, integral_pdf_time_trunc,
                                         integral_pdf_time_trunc_p,
                                         integral_pdf_time_trunc_c)



#%% ***************************************************************************


def clambdaj(theta, j, t, x, y, m, bk, ta, mc):
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
    return s



#%% ***************************************************************************


def clambdajGr(theta, j, t, x, y, m, bk, fv, dfv, ta, mc):
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
    
    part1_alpha = part1 * mc[:j] # m[:j]
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


def cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, ta, mc,
            voronoi=None):

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
        
    sk = A * np.exp(alpha * mc[j]) #m[j])
    return sk * gi * si



#%% ***************************************************************************


def cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv,
              ta, mc, voronoi=None):

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
    
    sk = A * np.exp(alpha * mc[j]) #m[j])
    fv = sk * gi * si
    dfv[0] = 0.
    dfv[1] = sk * gi  * si / A        * 2 * theta[1]
    dfv[2] = sk * gic * si            * 2 * theta[2]
    dfv[3] = sk * gi  * si * mc[j]    * 2 * theta[3] # m[j]
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
                  23.31388888888889,
                  23.378854166666665,
                  23.652256944444446])
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
    print(clambdaj(theta, j, t, x, y, m, bk, ta, 3.5))

    print('\nclambdajGr')
    fv = None
    dfv = [None]*8
    print(clambdajGr(theta, j, t, x, y, m, bk, fv, dfv, ta, 3.5))
    
    
    start_time = time.time()
    print('\ncintegj')
    print(cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, ta, 3.5))
    
    print('\ncintegjGr')
    fv = None
    dfv = [None]*8
    print(cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv, ta, 3.5))
    print(time.time()-start_time)

    