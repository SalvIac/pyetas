# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""


import numpy as np

from pyrisk.etas.etas8p.dist import dist, dist2, norm
import pyrisk.etas.etas8p.lambdaf0 as m0
import pyrisk.etas.etas8p.lambdaf1 as m1
import pyrisk.etas.etas8p.lambdaf2 as m2
import pyrisk.etas.etas8p.lambdaf3 as m3
import pyrisk.etas.etas8p.lambdaf4 as m4
import pyrisk.etas.etas8p.lambdaf5 as m5



#%% ***************************************************************************


def clinesearch(rdata, xOld, h, fv, verbose, ram, model_module):

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
    fv2 = cloglkhd(xNew, rdata, verbose, model_module)

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
            fv2 = cloglkhd(xNew, rdata, verbose, model_module)
    else:
        # stat30:
        ram3 = ram2*2.0
        for i in range(0, len(xOld)):
            xNew[i] = xOld[i] + ram3 * h[i]
        fv3 = cloglkhd(xNew, rdata, verbose, model_module)
        while not (fv3 > fv2):
            ram1 = ram2
            ram2 = ram3
            fv1 = fv2
            fv2 = fv3
            ram3 = ram2*2.0
            for i in range(0, len(xOld)):
                xNew[i] = xOld[i] + ram3 * h[i]
            fv3 = cloglkhd(xNew, rdata, verbose, model_module)

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
        fv = cloglkhd(xNew, rdata, verbose, model_module)
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
        fv = cloglkhd(xNew, rdata, verbose, model_module)
        if (fv2 < fv):
            ram = ram2
        return ram


#%% ***************************************************************************
# log-likelihood function of the model

def cloglkhd(tht, rdata, verbose, model_module):
    
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

    # extract polygon information
    px = rdata['rpoly']['px']
    py = rdata['rpoly']['py']
    npoly = len(py)

    # extract time period information
    tstart2 = rdata['tperiod']['study_start']
    tlength = rdata['tperiod']['study_end']

    # extract integral of spatial intensity over the obsevation window
    integ0 = rdata['integ0']

    fv1 = 0.
    fv2 = 0.

    for j in range(0, N):
        
        if fault[j] is None:
        
            if flag[j] == 1:
                s = model_module.clambdaj(tht, j, t, x, y, m, bk)
                if s > 1.0e-25:
                    fv1 += np.log(s)
                else:
                    fv1 -= 100.0
            fv2 += model_module.cintegj(tht, j, t, x, y, m, npoly, px, py, tstart2, tlength)
        
        else: ######################################################## fault
            num = fault[j]['xx'].shape[0] * fault[j]['xx'].shape[1]
            tht2 = tht.copy()
            tht2[1] = tht2[1]/num # A

            for jj, (xx, yy, dd) in enumerate(zip(fault[j]['xx'].flatten(),
                                                  fault[j]['yy'].flatten(),
                                                  fault[j]['depth'].flatten())):
                if flag[j] == 1:
                    
                    x2 = x.copy()
                    y2 = y.copy()
                    x2[j] = xx
                    y2[j] = yy
                    
                    s = model_module.clambdaj(tht2, j, t, x2, y2, m, bk)
                    if s > 1.0e-25:
                        fv1 += np.log(s)
                    else:
                        fv1 -= 100.0
        
                fv2 += model_module.cintegj(tht2, j, t, x2, y2, m, npoly, px, py, tstart2, tlength)

    fv2 += tht[0]**2 * (integ0)
    fv = -fv1 + fv2

    if verbose == 1:
        print("Function Value = {:8.5f}  {:7.2f}  {:7.2f}".format(fv, -fv1, fv2))

    return fv


#%% ***************************************************************************


def cloglkhdGr(tht, rdata, verbose, fv, dfv, model_module):
    
    print('................fault...............')
    
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

    # extract polygon information
    px = rdata['rpoly']['px']
    py = rdata['rpoly']['py']
    npoly = len(py)

    # extract time period information
    tstart2 = rdata['tperiod']['study_start']
    tlength = rdata['tperiod']['study_end']

    # extract integral of spatial intensity over the obsevation window
    integ0 = rdata['integ0']
    
    fv1 = 0.
    fv2 = 0.
    df1 = [0.]*len(tht)
    df2 = [0.]*len(tht)
    fv1temp = None
    g1temp = [None]*len(tht)
    fv2temp = None
    g2temp = [None]*len(tht)

    for j in range(0, N):
        
        if fault[j] is None:
                
            if flag[j] == 1:
                fv1temp, g1temp = model_module.clambdajGr(tht, j, t, x, y, m, bk, fv1temp, g1temp)
                if fv1temp > 1.0e-25:
                    fv1 += np.log(fv1temp)
                else:
                    fv1 -= 100.0
    
                for i in range(0, len(tht)):
                    df1[i] += g1temp[i] / fv1temp
    	
            fv2temp, g2temp = model_module.cintegjGr(tht, j, t, x, y, m, npoly, px, py, tstart2, tlength, fv2temp, g2temp)
            fv2 += fv2temp
            for i in range(0, len(tht)):
    	        df2[i] += g2temp[i]
        
        else: ################### fault
            num = fault[j]['xx'].shape[0] * fault[j]['xx'].shape[1]
            tht2 = tht.copy()
            tht2[1] = tht2[1]/num # A

            for jj, (xx, yy, dd) in enumerate(zip(fault[j]['xx'].flatten(),
                                                  fault[j]['yy'].flatten(),
                                                  fault[j]['depth'].flatten())):
                if flag[j] == 1:
                    
                    x2 = x.copy()
                    y2 = y.copy()
                    x2[j] = xx
                    y2[j] = yy
                    
                    fv1temp, g1temp = model_module.clambdajGr(tht2, j, t, x2, y2, m, bk, fv1temp, g1temp)
                    if fv1temp > 1.0e-25:
                        fv1 += np.log(fv1temp)
                    else:
                        fv1 -= 100.0
        
                    for i in range(0, len(tht)):
                        df1[i] += g1temp[i] / fv1temp
        	
                fv2temp, g2temp = model_module.cintegjGr(tht2, j, t, x2, y2, m, npoly, px, py, tstart2, tlength, fv2temp, g2temp)
                fv2 += fv2temp
                for i in range(0, len(tht2)):
        	        df2[i] += g2temp[i]
                
                
    fv2 += tht[0]**2 * integ0
    df2[0] = integ0 * tht[0] * 2

    fv = -fv1 + fv2
    
    for i in range(0, len(tht)):
        dfv[i] = -df1[i] + df2[i]
        
    if verbose == 1:
        print("Function Value = {:8.2f}  {:7.2f}  {:7.2f}".format(fv, -fv1, fv2))
        for i in range(0, len(tht)):
            print("Gradiant [{:d}] = {:8.2f}    theta[{:d}] = {:2.8f}".format(i + 1, dfv[i], i + 1, tht[i]))
    
    return fv, dfv


# *******************************************************************************


def cfit(theta, rdata, ihess, rverbose, model_module):

    # SEXP estimate, fvout, dfvout, aic, hess, out
    
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
    
    fv, g = cloglkhdGr(tht, rdata, verbose, fv, g, model_module)

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
                        print("theta[{:d}] = {:2.8f}\t gradient[{:d}] = {:8.4f}\n".format(
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
            ramda = clinesearch(rdata, tht, s, ed, verbose, ramda, model_module)
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
            fv, g = cloglkhdGr(tht, rdata, verbose, fv, g, model_module)

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
                        print("theta[{:d}] = {:2.8f}\t gradient[{:d}] = {:8.4f}\n".format(
                              i + 1, tht[i]**2, i + 1, g[i]))
                out = dict(estimate = estimP,
                           fvout = fvP,
                           dfvout = dfvP,
                           aic = aicP,
                           hess = hessP)
                return out
    return 0


#%% ***************************************************************************


def etasfit(theta, revents, rpoly, tperiod, integ0, ihess, verbose, ndiv, eps,
            cxxcode, nthreads, model):
    tht = {j[0]: np.sqrt(j[1]) for j in theta.items()}
    # if False: #cxxcode
    #     cfit_res = cxxfit(tht, revents, rpoly, tperiod, integ0, ihess,
    #                   int(ndiv), eps, bool(verbose), int(nthreads))
    
    rdata = dict(revents=revents, rpoly=rpoly, tperiod=tperiod, integ0=integ0)
    if model == 0: # import functions for model 0
        print('imported model 0 - space-independent')
        model_module = m0
    elif model == 1: # import functions for model 1
        model_module = m1
    elif model == 2: # import functions for model 2
        model_module = m2
    elif model == 3: # import functions for model 3
        model_module = m3
    elif model == 4: # import functions for model 4
        model_module = m4
    elif model == 5: # import functions for model 5
        model_module = m5
    else:
        raise Exception('invalid model')
    
    cfit_res = cfit(tht, rdata, ihess, int(verbose), model_module)
    
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
         
