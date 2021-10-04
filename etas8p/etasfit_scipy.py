# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""


import numpy as np
from scipy import optimize # minimize, LinearConstraint, NonlinearConstraint, Bounds

from pyrisk.etas.etas8p.etasfit import cloglkhd, cloglkhdGr
# import pyrisk.etas.etas8p.lambdaf0 as m0
import pyrisk.etas.etas8p.lambdaf1 as m1
import pyrisk.etas.etas8p.lambdaf2 as m2
import pyrisk.etas.etas8p.lambdaf3 as m3
import pyrisk.etas.etas8p.lambdaf4 as m4
import pyrisk.etas.etas8p.lambdaf5 as m5
import pyrisk.etas.etas8p.lambdaf6 as m6
import pyrisk.etas.etas8p.lambdaf6f as m6f
import pyrisk.etas.etas8p.lambdaf7 as m7
import pyrisk.etas.etas8p.lambdaf7f as m7f



#%% ***************************************************************************


def grad(tht, rdata, verbose, model_module, voronoi):
    fv, g = cloglkhdGr(tht, rdata, verbose, None, [None]*len(tht),
                       model_module, voronoi)
    return np.array(g)


def etasfit_scipy(theta, revents, rpoly, tperiod, integ0, m0, ihess, verbose,
                  ndiv, eps, model, impbounds=None, modified_etas=False,
                  voronoi=None):
    rdata = dict(revents=revents, rpoly=rpoly, tperiod=tperiod, integ0=integ0, mmin=m0)
    x0 = [np.sqrt(j[1]) for j in theta.items()]
    linear_constraint = []
    
    # # fix alpha=beta
    # A = np.zeros((8,8))
    # A[3,3] = 1.
    # lb = ub = np.array([0.,0.,0.,x0[3],0.,0.,0.,0.])
    # linear_constraint = [optimize.LinearConstraint(A, lb, ub)]

    # # fix q=1.5
    # A = np.zeros((8,8))
    # A[6,6] = 1.
    # lb = ub = np.array([0.,0.,0.,0.,0.,0.,x0[6],0.])
    # linear_constraint = [optimize.LinearConstraint(A, lb, ub)]
    
    # # fix q=1.5 and alpha=beta
    # A = np.zeros((8,8))
    # A[3,3] = 1.
    # A[6,6] = 1.
    # lb = ub = np.array([0.,0.,0.,x0[3],0.,0.,x0[6],0.])
    # linear_constraint = [optimize.LinearConstraint(A, lb, ub)]
    
    fault_mode = not np.all([np.array(revents["fault"]) == None])
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
        
    else:
        raise Exception('invalid model')
    
    # default constraints
    bounds_seq = [(lb, ub) for lb, ub in zip(
                                            [-np.inf]*len(list(theta.keys())),
                                            [+np.inf]*len(list(theta.keys())))]
    print(bounds_seq)
    
    if impbounds is not None: # update the constraints in case they are specified
        if not isinstance(impbounds, dict): raise Exception("impbounds must be a dict!")
    
        for key in impbounds.keys():
            idx = list(theta.keys()).index(key)
            if isinstance(impbounds[key], list): # for ranges
                if impbounds[key][0] >= impbounds[key][1]:
                    raise Exception("if you want to specify a range for "+key+
                                    " the first element of the list must be < the second element of the list")
                bounds_seq[idx] = (np.sqrt(impbounds[key][0]), np.sqrt(impbounds[key][1]))
            else: # fix specific value (more or less)
                if impbounds[key] == 0:
                    bounds_seq[idx] = (np.sqrt(impbounds[key]),
                                       np.sqrt(impbounds[key]+1e-6))
                else:
                    bounds_seq[idx] = (np.sqrt(impbounds[key]-1e-6),
                                       np.sqrt(impbounds[key]+1e-6))            
        print(bounds_seq)


    # if list(theta.keys()) == ['mu', 'A', 'c', 'alpha', 'p', 'D', 'q', 'gamma']: # theta['alpha']-1e-6
    #     # bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip(
    #     #     [  0.01, 0.0001,   1e-9,   0.01, 0.0001,   1e-9, 1.0001,   0.01],
    #     #     [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])]
    #     # bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip(
    #     #     [  0.01, 0.01,   1e-9, theta['alpha']-1e-6, 0.5,   1e-7, 1.0001, 0.0001],
    #     #     [    2.,   2.,    0.1, theta['alpha']+1e-6,  2.,     1., np.inf, np.inf])]
    #     # bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip(
    #     #     [  0.01, 0.01,   1e-9, 0.01, 0.5,   1e-7, theta['q']-1e-6, 0.0001],
    #     #     [    2.,   2.,    0.1,   3.,  2.,     1., theta['q']+1e-6, np.inf])]
    #     bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip(
    #         [  0.01, 0.01,   1e-9, theta['alpha']-1e-6, 0.5,   1e-7, theta['q']-1e-6, 0.0001],
    #         [    2.,   2.,    0.1, theta['alpha']+1e-6,  2.,     1., theta['q']+1e-6, np.inf])]
    #     # bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip(
    #     #    [theta['mu']-1e-6, 0.0001,   1e-9,   0.01, 0.0001,   1e-9, 1.0001,   0.01],
    #     #    [theta['mu']+1e-6, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])]
    #     # bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip(
    #     #     [theta['mu']-1e-6, 0.0001,   1e-9, theta['alpha']-1e-6, 0.0001,   1e-9, 1.0001,   0.01],
    #     #     [theta['mu']+1e-6, np.inf, np.inf, theta['alpha']+1e-6, np.inf, np.inf, np.inf, np.inf])]

    # elif list(theta.keys()) == ['mu', 'A', 'c', 'alpha', 'p', 'D', 'q']:
    #     bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip([  0.01,   0.01,   1e-9,   0.01, 1.0001,   1e-9, 1.0001],
    #                                                               [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])]
    # elif list(theta.keys()) == ['mu', 'A', 'c', 'alpha', 'p', 'D']:
    #     bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip([  0.01,   0.01,   1e-9,   0.01, 1.0001,   1e-9],
    #                                                               [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])]

    # elif list(theta.keys()) == ['mu', 'A', 'c', 'alpha', 'p']:
    #     bounds_seq = [(np.sqrt(lb),np.sqrt(ub)) for lb, ub in zip([  0.01,   0.01,   1e-9,   0.01, 1.0001],
    #                                                               [np.inf, np.inf, np.inf, np.inf, np.inf])]

    bounds = optimize.Bounds([b[0] for b in bounds_seq], [b[1] for b in bounds_seq])

    # res = optimize.minimize(cloglkhd, x0,
    #                         args=(rdata, verbose, model_module, voronoi), 
    #                         method='trust-constr', jac=grad,
    #                         hess=optimize.BFGS(), bounds=bounds,
    #                         constraints=linear_constraint,
    #                         options={'xtol': 1e-06, 'gtol': 1e-04,
    #                                   'barrier_tol': 1e-05, 'verbose': 1})
    res = optimize.minimize(cloglkhd, x0,
                            args=(rdata, verbose, model_module, voronoi), 
                            method='BFGS', jac=grad, #bounds=bounds,
                            options={'disp': True}) # , "maxcor": int(1e3)
    
    # res = optimize.minimize(cloglkhd, x0,
    #                         args=(rdata, verbose, model_module, voronoi),
    #                         jac=grad, bounds=bounds,
    #                         options={'ftol': 0.2e-11, 'gtol': 1e-6})
    
    # # global
    # results_DE = optimize.differential_evolution(cloglkhd, bounds_seq, args=(rdata, verbose))
    # results_shgo = optimize.shgo(cloglkhd, bounds_seq, args=(rdata, verbose))
    # results_DA = optimize.dual_annealing(cloglkhd, bounds_seq, x0=x0, args=(rdata, verbose))
    
    # unconstrained
    # res = optimize.minimize(cloglkhd, x0, args=(rdata, verbose, model_module, voronoi),
    #                         method='nelder-mead',
    #                         options={'xatol': 1e-8, 'disp': True})
    # res = optimize.minimize(cloglkhd, x0, args=(rdata, verbose), method='L-BFGS-B', jac=grad, options={'disp': True})
    # res = optimize.minimize(cloglkhd, x0, args=(rdata, verbose), method='BFGS', jac=grad, options={'disp': True})
    # res = optimize.minimize(cloglkhd, x0, args=(rdata, verbose), method='Newton-CG', jac=grad, options={'disp': True})
    # res = optimize.minimize(cloglkhd, x0, args=(rdata, verbose), method='trust-ncg', jac=grad, options={'disp': True})

    
    H = np.zeros((len(x0),len(x0)))
    tht = np.array(res.x)
    avcov = (1/4 * np.matmul(np.matmul(np.diag(1/tht), H), np.diag(1/tht))).tolist()

    cfit_res = dict()
    cfit_res['estimate'] = [t**2 for t in tht]
    cfit_res['avcov'] = avcov
    cfit_res['loglik'] = -res.fun
    cfit_res['aic'] = 2 * (-cfit_res['loglik'] + len(tht))
    try:
        cfit_res['gradient'] = res.grad
    except:
        cfit_res['gradient'] = np.zeros(len(x0))
    cfit_res['ihessian'] = np.zeros((len(x0),len(x0)))
    cfit_res['res'] = res
    return cfit_res
         

#%%


if __name__ == "__main__":

    import time
    import sys
    sys.path.append('C:\\Users\\Salvatore\\Dropbox\\SalvIac')
    from pyrisk.etas.etas8p.voronoi import get_voronoi
    from myutils.utils_pickle import load_pickle, save_pickle
    
    theta = load_pickle('test/param1')
    rdata = load_pickle('test/rdata')
    revents = rdata['revents']
    rpoly = rdata['rpoly']
    tperiod = rdata['tperiod']
    integ0 = rdata['integ0']
    ihess = load_pickle('test/ihess')
    rverbose = verbose = 1

    ndiv = 1000
    eps=1e-06
    # cxxcode=False
    # nthreads=1  
    m0 = 3.5
    revents["fault"] = [None]*len(revents["mm"])
    rdata['mmin'] = m0
    tht = [np.sqrt(j[1]) for j in theta.items()]
    
    # # voronoi tassellation
    # time1 = time.time()
    # points, areas, _ = get_voronoi(rpoly, min_prec=0.05)
    # voronoi = {"points": points,
    #            "areas": areas}
    # print("create vonoroi", time.time()-time1)




    # ##########################################################################
    # print('\netasfit constrained with scipy minimize model 6')
    # time1 = time.time()
    # print(theta)
    # res = etasfit_scipy(theta, revents, rpoly, tperiod, integ0, m0, ihess,
    #                     verbose, ndiv, eps, model=6, impbounds={"q":1.5},
    #                     modified_etas=False, voronoi=None)
    # print(time.time()-time1)
    # print(res)
    
    


    ##########################################################################
    print('\netasfit with scipy minimize model 6')
    time1 = time.time()
    print(theta)
    res = etasfit_scipy(theta, revents, rpoly, tperiod, integ0, m0, ihess,
                        verbose, ndiv, eps, model=6, modified_etas=False,
                        voronoi=None)
    print(time.time()-time1)
    print(res)
    



    ##########################################################################
    # print('\netasfit with scipy minimize model 5')
    # time1 = time.time()
    # print(theta)
    # res = etasfit(theta, revents, rpoly, tperiod, integ0,
    #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=5)
    # print(time.time()-time1)
    # print(res)
    # results = [0.555266260648366, 0.1864583473888933, 0.047195657755165474,
                # 2.7050217320762324, 1.1547553549260672, 0.015949991076632867,
                # 2.317120770715824, 0.023512565102052223]




    ##########################################################################
    # print('\netasfit with scipy minimize model 4')
    # time1 = time.time()
    # del theta['gamma']
    # ihess = ihess[:7,:7]
    # print(theta)
    # res = etasfit(theta, revents, rpoly, tperiod, integ0,
    #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=4)
    # print(res)
    # print(time.time()-time1)




    ##########################################################################
    # print('\netasfit with scipy minimize model 3')
    # time1 = time.time()
    # del theta['gamma']
    # ihess = ihess[:7,:7]
    # print(theta)
    # res = etasfit(theta, revents, rpoly, tperiod, integ0,
    #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=3)
    # print(res)
    # print(time.time()-time1)




    ##########################################################################
    # print('\netasfit with scipy minimize model 2')
    # time1 = time.time()
    # del theta['gamma']
    # del theta['q']
    # ihess = ihess[:6,:6]
    # print(theta)
    # res = etasfit(theta, revents, rpoly, tperiod, integ0,
    #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=2)
    # print(res)
    # print(time.time()-time1)




    ##########################################################################
    # print('\netasfit with scipy minimize model 1')
    # time1 = time.time()
    # del theta['gamma']
    # del theta['q']
    # ihess = ihess[:6,:6]
    # print(theta)
    # res = etasfit(theta, revents, rpoly, tperiod, integ0,
    #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=1)
    # print(res)
    # print(time.time()-time1)




    ##########################################################################
    # print('\netasfit with scipy minimize model 0')
    # time1 = time.time()
    # del theta['gamma']
    # del theta['q']
    # del theta['D']
    # ihess = ihess[:5,:5]
    # print(theta)
    # res = etasfit(theta, revents, rpoly, tperiod, integ0,
    #               ihess, verbose, ndiv, eps, cxxcode, nthreads, model=0)
    # print(res)
    # print(time.time()-time1)
    



