# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""

import time
import warnings
from tqdm import tqdm
import numpy as np
import scipy.spatial
import pandas as pd
import datetime

from pyrisk.etas.etas8p.catalog import CatalogueEtas
from pyrisk.etas.etas8p.decluster import decluster
from pyrisk.etas.etas8p.etasfit import etasfit
from pyrisk.etas.etas8p.dist import dist
from myutils.utils_pickle import load_pickle, save_pickle



class EtasMle():
        
    def __init__(self, catalog, param0=None, bwdf=None, pb=None, integ0=None,
                 verbose=True, plot_it=False, ndiv=1000,
                 rel_tol=1e-03, eps=1e-06, cxxcode=True, nthreads=1,
                 model=5):
        
        ptm = time.time()
        revents = catalog.revents
        rpoly = catalog.rpoly
        rtperiod = catalog.rtperiod
        m0 = catalog.mag_threshold
        win = catalog.region_win
        
        # initial prameter values
        if param0 is None:
            mu0 = revents.shape[0]/(4*(rtperiod['study_end']-rtperiod['study_start']) * win.area) # N/(4*T*S)
            param0 = dict(mu=mu0, A=0.01, c=0.01, alpha=1., p=1.3, D=0.01, q=2., gamma=1.)
            if catalog.dist_unit == "km":
                param0["D"] = 111**2 * param0["D"]
            if verbose:
                print("using non-informative initial parameter values:\n")
                print(param0)
            print("the algorithm is very sensitive to the choice of starting point")


        # bandwidths for smoothness and integration
        if bwdf is None:
            bwdf = 0.05
        if catalog.dist_unit == "km":
            bwdf = 6371.3 * np.pi / 180 * bwdf
        rbwd = bwdf*np.ones(revents.shape[0])
        
        # check initial values for the model parameters
        if any([value < 0. for _, value in param0.items()]): # len(list(param0)) != 8 or
            raise Exception("param0 must be a numeric vector of length 8 with positive components")

        param1 = param0.copy()
        thetar = np.empty((1,len(param0.keys()),))
        thetar[:] = np.nan
        asd = np.empty((1,len(param0.keys()),))
        asd[:] = np.nan
        par_names = list(param1.keys())
        loglikfv = np.zeros(revents.shape[0])
        loglikfv[:] = np.nan
        ihess = np.diag(np.full(len(param0.keys()),1.))
        events = {j[0]: j[1].to_list() for j in revents.iteritems()}
    
    
        tht = {j[0]: np.sqrt(j[1]) for j in theta.items()}
        cbkg = cdeclust(tht, rbwd, revents, rpoly, rtperiod, model)
        events = bkg['revents']
        integ0 = bkg['integ0']
        bk = np.array(events['bkgd'])
        pb = np.array(events['prob'])
    
        
    
    
    
    
        if verbose:
            print("======================================================")
            print("\nintegral of background seismicity rate: ", integ0)
            print("======================================================")

        print("estimating:")
        opt = etasfit(param1, events, rpoly, rtperiod, integ0, ihess,
                      verbose, ndiv, eps, cxxcode, nthreads, model)
        thetar[0,:] = opt['estimate']
        loglikfv[0] = opt['loglik']
        asd[0,:] = np.sqrt(np.diag(opt['avcov']))
        ihess = np.array(opt['ihessian'])
        
        if model==1:
            param1 = dict(mu=thetar[0,0], A=thetar[0,1], c=thetar[0,2],
                          alpha=thetar[0,3], p=thetar[0,4], D=thetar[0,5])
        elif model==2:
            param1 = dict(mu=thetar[0,0], A=thetar[0,1], c=thetar[0,2],
                          alpha=thetar[0,3], p=thetar[0,4], D=thetar[0,5])
        elif model==3:
            param1 = dict(mu=thetar[0,0], A=thetar[0,1], c=thetar[0,2],
                          alpha=thetar[0,3], p=thetar[0,4], D=thetar[0,5],
                          q=thetar[0,6])
        elif model==4:
            param1 = dict(mu=thetar[0,0], A=thetar[0,1], c=thetar[0,2],
                          alpha=thetar[0,3], p=thetar[0,4], D=thetar[0,5],
                          q=thetar[0,6])
        elif model==5:
            param1 = dict(mu=thetar[0,0], A=thetar[0,1], c=thetar[0,2],
                          alpha=thetar[0,3], p=thetar[0,4], D=thetar[0,5],
                          q=thetar[0,6], gamma=thetar[0,7])
            
        if verbose:
            print("======================================================")
            print("MLE:\n")
            print(param1)
            print("======================================================")
            # other formulation of etas have the productivity k0 instead of A
            # k0 = 1/np.pi*A*(p-1)*c**(p-1)*(q-1)*D**(q-1)
        

        if verbose:
            print("Execution time:", time.time() - ptm,'s')
        
        # self.catalog.revents = pd.DataFrame(events)
        self.catalog = catalog
        self.param1 = param1
        self.revents = events
        self.param = param1
        self.bk = bk
        self.pb = pb
        self.opt = opt
        self.bwd = bwd
        self.thetar = thetar
        self.loglikfv = loglikfv
        self.asd = asd
        self.integ0 = integ0
        self.ndiv = ndiv
        self.itr = None
        self.exectime = time.time()-ptm
    
    
    
    def print(self):
    
        flag = np.array(self.revents['flag'])
        mm = np.array(self.revents['mm'])[flag == 1]
        bt = 1. / mm.mean()
        asd_bt = bt**2 / mm.shape[0]
        param = self.param
        std = {key: self.asd[self.itr,i]  for i, key in enumerate(param.keys())}
        
        param1 = pd.DataFrame(dict(beta=[bt, asd_bt]), index=['Estimate', 'StdErr'])
        param2 = pd.DataFrame(self.param, index=['Estimate'])
        param3 = pd.DataFrame(std, index=['StdErr'])
        param23 = pd.concat([param2, param3], axis=0)
        params = pd.concat([param1, param23], axis=1)
        
        string = "ETAS model: fitted using iterative stochastic declustering method\n" + \
                 "converged after " + str(self.itr+1) + " iterations: elapsed execution time " + \
                 str(round(self.exectime/60, 2)) + " minutes\n\n" + \
                 "ML estimates of model parameters:\n" + str(params) + \
                 "\nDeclustering probabilities:\n" + str(pd.DataFrame(self.pb).describe().transpose()) + \
                 "\nlog-likelihood: " + str(self.opt['loglik']) + \
                 "\tAIC: " + str(self.opt['aic']) + "\n"
        return dict(string=string, params=params, loglik=self.opt['loglik'], aic=self.opt['aic'])


    def plot(self, which="est"):
        
        import matplotlib.pyplot as plt
        
        if which == "loglik":
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(np.arange(1,self.itr+1), self.loglikfv[0:self.itr], marker='o')
            ax.set_title("log-likelihood function of the model")
            plt.show()
        
        elif which == "est":
            keys = list(self.param.keys())
            data = self.thetar[0:self.itr+1,:]
            plt.figure()
            fig, axs = plt.subplots(4, 2)
            for i in range(0,8):
                axs[int(i/2)][i%2].set_title(keys[i])
                axs[int(i/2)][i%2].plot(np.arange(1,self.itr+2), data[:,i], marker='o')
            fig.suptitle('estimates of the model parameters')
            
        elif which == "rates":
            print('      rates.inter(x$param, x$object, x$bwd) ')
        else:
            raise Exception("Wrong type")
            




def backgr(theta, rbwd, revents, rpoly, tperiod, model):

    # extract events
    t = revents['tt']
    x = revents['xx']
    y = revents['yy']
    m = revents['mm']
    bk = revents['bkgd']
    pb = np.array(revents['prob'])
    lam = revents['lambd']
    t_arr = np.array(t)
    x_arr = np.array(x)
    y_arr = np.array(y)
    m_arr = np.array(m)
    # events = [revents[j] for j in revents.keys()]
    N = len(revents['tt'])

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
        r0 = dist(x[i], y[i], x_arr, y_arr)
        dg = dGauss(np.array(r0), bwd)
        s = np.sum(pb * np.array(dg))
        bk[i] = s / (tlength - tstart2)
    revents['bkgd'] = bk

    temp = [None]*N
    for i in range(0, N):
        temp[i] = polyinteg(pGauss, np.array([bwd[i]]), npoly, px, py, x[i], y[i])
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
        else:
            raise Exception('problem with model module')
            
    s = np.sum(pb * np.array(temp))
    revents['prob'] = pb
    revents['lambd'] = lam

    out = dict(revents=revents, integ0=s)
    return out    





