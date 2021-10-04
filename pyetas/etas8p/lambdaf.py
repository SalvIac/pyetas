# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""


import numpy as np

from pyrisk.etas.etas8p.dist import dist2
from pyrisk.etas.etas8p.poly import polyinteg, fr, dgamma_fr, dD_fr, dq_fr


#%% ***************************************************************************

# def fr(r, w):
#     gamma = w[0]
#     D = w[1]
#     q = w[2]
#     mag = w[3]
#     sig = D * np.exp(gamma * mag)
#     return (1 - (1 + r**2 / sig)**(1 - q)) / (2 * np.pi)

# def dgamma_fr(r, w):
#     gamma = w[0]
#     D = w[1]
#     q = w[2]
#     mag = w[3]
#     sig = D * np.exp(gamma * mag)
#     return (1 - q) * (1 + r**2 / sig)**(-q) * mag * r**2 / sig / (2 * np.pi)

# def dD_fr(r, w):
#     gamma = w[0]
#     D = w[1]
#     q = w[2]
#     mag = w[3]
#     sig = D * np.exp(gamma * mag)
#     return (1 - q) * (1 + r**2 / sig)**(-q) / D * r**2 / sig / (2 * np.pi)

# def dq_fr(r, w):
#     gamma = w[0]
#     D = w[1]
#     q = w[2]
#     mag = w[3]
#     sig = D*np.exp( gamma*mag )
#     return (1 + r**2 / sig)**(1 - q) * np.log(1 + r**2 / sig) / (2 * np.pi)


#%% ***************************************************************************


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


#%% ***************************************************************************


def clambdajGr(theta, j, t, x, y, m, bk, fv, dfv):
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
    part1 = np.exp(alpha * m[:j])
    
    delta = t[j] - t[:j]
    part2 = (p - 1)/c * np.power(1 + delta / c, - p)
    
    sig = D * np.exp(gamma * m[:j])
    r2 = dist2(x[j], y[j], x[:j], y[:j])
    part3 = (q - 1)/(sig * np.pi) * np.power(1 + r2/sig, - q)
    
    part_s = A * part1 * part2 * part3
    s = mu * bk[j] + np.sum(part_s)

    sg1 = bk[j]
    
    sg2 = np.sum(part1 * part2 * part3)
        
    part2_c = part2 * (-1/c - p/(c + delta) + p/c)
    sg3 = A * np.sum(part1 * part2_c * part3)
    
    part1_alpha = part1 * m[:j]
    sg4 = A * np.sum(part1_alpha * part2 * part3)
    
    part2_p = part2 * (1/(p - 1) - np.log(1 + delta/c))
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



def clambdajGr_old(theta, j, t, x, y, m, bk, fv, dfv):
    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    
    # double part1, part2, part3, part1_alpha, part2_c, part2_p, part3_d, part3_q,
    # part3_gamma, delta, sig, r2, sg1,
    sg1 = bk[j]
    sg2 = 0.
    sg3 = 0.
    sg4 = 0.
    sg5 = 0.
    sg6 = 0.
    sg7 = 0.
    sg8 = 0.

    s = mu * bk[j]

    for i in range(0, j):
        part1 = np.exp(alpha * m[i])
        
        delta = t[j] - t[i]
        part2 = (p - 1)/c * np.power(1 + delta / c, - p)
        
        sig   = D * np.exp(gamma * m[i])
        r2 = dist2(x[j], y[j], x[i], y[i])
        part3 = (q - 1)/(sig * np.pi) * np.power(1 + r2/sig, - q)
        
        s += A * part1 * part2 * part3
        sg2 += part1 * part2 * part3
        
        part2_c = part2 * (-1/c - p/(c + delta) + p/c)
        sg3 += A * part1 * part2_c * part3
        
        part1_alpha = part1 * m[i]
        sg4 += A * part1_alpha * part2 * part3
        
        part2_p = part2 * (1/(p - 1) - np.log(1 + delta/c))
        sg5 += A * part1 * part2_p * part3
        
        part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)))
        sg6 += A * part1 * part2 * part3_d
        
        part3_q = part3 * (1/(q - 1) - np.log(1 + r2/sig))
        sg7 += A * part1 * part2 * part3_q
        
        part3_gamma = part3 * (-m[i] + q * m[i] * (1 - 1/(1 + r2/sig)))
        sg8 += A * part1 * part2 * part3_gamma
          
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


def cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength):

    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    
    # double ttemp, ttemp1, ttemp2, gi, gi1, gi2, w[4], si, sk
    w = [None]*4
    
    if t[j] > tstart2:
        ttemp = tlength - t[j]
        gi  = 1 - (1 + ttemp/c)**(1 - p)
    else:
        ttemp1 = tstart2 - t[j]
        ttemp2 = tlength - t[j]
        
        gi1  = 1 - (1 + ttemp1/c)**(1 - p)
        gi2  = 1 - (1 + ttemp2/c)**(1 - p)
        gi  = gi2 - gi1

    w = np.array([ gamma, D, q, m[j] ])
    si = polyinteg(fr, w, npoly, px, py, x[j], y[j])
    sk = A * np.exp(alpha * m[j])
    return sk * gi * si


#%% ***************************************************************************


def cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv):

    # extract model parameters
    mu = theta[0]**2
    A = theta[1]**2
    c = theta[2]**2
    alpha = theta[3]**2
    p = theta[4]**2
    D = theta[5]**2
    q = theta[6]**2
    gamma = theta[7]**2
    
    # double ttemp, ttemp1, ttemp2, gi, gi1, gi2, gic, gic1, gic2, gip, gip1,
    # gip2, w[4], si, sid, siq, sigamma, sk
    w = [None]*4
    
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

    w = np.array([ gamma, D, q, m[j] ])
    si      = polyinteg(fr, w, npoly, px, py, x[j], y[j])
    sid     = polyinteg(dD_fr, w, npoly, px, py, x[j], y[j])
    siq     = polyinteg(dq_fr, w, npoly, px, py, x[j], y[j])
    sigamma = polyinteg(dgamma_fr, w, npoly, px, py, x[j], y[j])
    
    sk = A * np.exp(alpha * m[j])
    fv      = sk * gi * si
    dfv[0] = 0
    dfv[1] = sk * gi  * si / A        * 2 * theta[1]
    dfv[2] = sk * gic * si            * 2 * theta[2]
    dfv[3] = sk * gi  * si * m[j]     * 2 * theta[3]
    dfv[4] = sk * gip * si            * 2 * theta[4]
    dfv[5] = sk * gi  * sid           * 2 * theta[5]
    dfv[6] = sk * gi  * siq           * 2 * theta[6]
    dfv[7] = sk * gi  * sigamma       * 2 * theta[7]
    return fv, dfv



#%% ***************************************************************************


# SEXP clambdax(SEXP rt,
#              SEXP rx,
#              SEXP ry,
#              SEXP theta,
#              SEXP revents)
# {
#   SEXP dim

#   # extract events
#   PROTECT(dim = allocVector(INTSXP, 2))
#   dim = getAttrib(revents, R_DimSymbol)
#   int N = INTEGER(dim)[0]
#   double *events = REAL(revents)
#   double t[N], x[N], y[N], m[N]
#   for (int i = 0 i < N i++)
#     {
#       t[i] = events[i]
#       x[i] = events[N + i]
#       y[i] = events[2 * N + i]
#       m[i] = events[3 * N + i]
#     }

#   # extract model parameters
#   double *tht = REAL(theta)
#   double #mu = tht[0] * tht[0],
#     A     = tht[1] * tht[1],
#     c     = tht[2] * tht[2],
#     alpha = tht[3] * tht[3],
#     p     = tht[4] * tht[4],
#     D     = tht[5] * tht[5],
#     q     = tht[6] * tht[6],
#     gamma = tht[7] * tht[7]

#   # extract arguments
#   double tt = *REAL(rt), xx = *REAL(rx), yy = *REAL(ry)

#   double part1, part2, part3, delta, sig, r2

#   double s = 0

#   int i = 0
#   while (t[i] < tt && i < N)
#     {
#       part1 = exp(alpha * m[i])

#       delta = tt - t[i]
#       part2 = (p - 1)/c * pow(1 + delta/c, - p)

#       sig = D * exp(gamma * m[i])
#       r2 = dist2(xx, yy, x[i], y[i])
#       part3 = (q - 1) / (sig * PI) * pow(1 + r2 / sig, - q)

#       s += A * part1 * part2 * part3
#       i++
#     }

#   SEXP out
#   PROTECT(out = allocVector(REALSXP, 1))
#   double *outP = REAL(out)
#   *outP = s
#   UNPROTECT(2)
#   return out
# }


#%% ***************************************************************************


if __name__ == "__main__":
    
    theta = [np.sqrt(j[1]) for j in {'mu': 0.11459513439367879,
                                    'A': 0.01,
                                    'c': 0.01,
                                    'alpha': 1,
                                    'p': 1.3,
                                    'D': 0.01,
                                    'q': 2,
                                    'gamma': 1}.items()]
    j = 5
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
                 1.0,
                 0.20000000000000018,
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
    
    import time
    
    # print('\nclambdaj')
    # # test ok with c++ 0.00594617 
    # print(clambdaj(theta, j, t, x, y, m, bk))
    # time1 = time.time()
    # for i in range(10000):
    #     clambdaj(theta, j, t, x, y, m, bk)
    # print(time.time()-time1)
    
    # print('\nclambdajGr')
    # fv = None
    # dfv = [None]*8
    # # test ok with c++
    # # 0.0351299 2.10445e-06 6.24062e-07 1.87068e-07 -6.24517e-07 2.08408e-06 -1.2921e-06 1.8524e-07 0.00594617        
    # print(clambdajGr(theta, j, t, x, y, m, bk, fv, dfv))
    # time1 = time.time()
    # for i in range(10000):
    #     clambdajGr(theta, j, t, x, y, m, bk, fv, dfv)
    # print(time.time()-time1)
    # print('clambdajGr_old')
    # print(clambdajGr_old(theta, j, t, x, y, m, bk, fv, dfv))
    # time1 = time.time()
    # for i in range(10000):
    #     clambdajGr_old(theta, j, t, x, y, m, bk, fv, dfv)
    # print(time.time()-time1)
    
    print('\ncintegj')
    # # test ok with c++ 0.0237023
    print(cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength))
    # time1 = time.time()
    # for i in range(10000):
    #     cintegj(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength)
    # print(time.time()-time1)
    
    print('\ncintegjGr')
    fv = None
    dfv = [None]*8
    # # test ok with c++
    # # 0 0.474047 -0.020691 0.0474047 0.0541212 -0.000495219 0.000454441 -4.95219e-05 
    print(cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv))
    # time1 = time.time()
    # for i in range(10000):
    #     cintegjGr(theta, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv, dfv)
    # print(time.time()-time1)



