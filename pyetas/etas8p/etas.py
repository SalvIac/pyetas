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
the core of the code is similar to:
https://github.com/jalilian/ETAS/tree/master/R
http://bemlar.ism.ac.jp/zhuang/software.html
with the changes explained in Iacoletti et al. 2021
"""

import time
import warnings
from tqdm import tqdm
import numpy as np
import scipy.spatial
import pandas as pd
import datetime
from copy import deepcopy

from openquake.hmtk.seismicity.occurrence.aki_maximum_likelihood import AkiMaxLikelihood

from pyetas.etas8p.catalog import CatalogueEtas
from pyetas.etas8p.decluster import decluster
from pyetas.etas8p.etasfit import etasfit
from pyetas.etas8p.etasfit_scipy import etasfit_scipy # q=1.5 alpha=beta
from pyetas.etas8p.voronoi import get_voronoi
from myutils.utils_pickle import load_pickle, save_pickle



class Etas():
        
    def __init__(self, catalog, param0=None, bwd=None, nnp=5, bwm=0.05,
                 verbose=True, pb_fix=False, ndiv=1000, no_itr=11,
                 rel_tol=1e-03, eps=1e-06, use_scipy=False, model=6, fault=True,
                 impbounds=None, voronoi_mode=False, modified_etas=False):
        
        '''
        The function etas fits the ETAS model to a catalog of earthquakes. The function takes
        an object of class ‘catalog’ and calls the internal functions decluster and etasfit for
        simultaneous estimation of the background seismicity rate u(x, y) and the parameter vector
        using the iterative approach discussed in Section 3. An initial guess for the parameter vector
        = (μ, A, , c, p,D, , q) can be provided by the param0 argument. Thus, param0 needs to
        be a numeric vector of length 8. If param0 is not specified, as suggested by Ogata (1998),
        the default values μ = N/(4T|S|), A = 0.01, c = 0.01, = 1, p = 1.3, D = 0.01, q = 2 and
         = 1 are used as components of param0. Note that this is a crude initial estimate and does
        not guarantee the convergence of Algorithm 2. The algorithm is relatively sensitive to the
        choice of the initial estimate of and, if possible, providing more informative initial estimate
        is highly recommended.
        '''
        
        ptm = time.time()
        revents = catalog.revents
        rpoly = catalog.rpoly
        rtperiod = catalog.rtperiod
        m0 = catalog.mag_threshold
        win = catalog.region_win
        
        # b value
        cat = deepcopy(catalog)
        b, bstd = self._get_b_value(cat)
        print("b-value", b)
        
        if not fault:
            catalog.revents["fault"] = None

        voronoi = None
        if voronoi_mode:
            points, areas, _ = get_voronoi(rpoly, min_prec=0.05)
            voronoi = {"points": points,
                       "areas": areas}
        
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
        if bwd is None:
            if catalog.dist_unit == "km":
                bwm = 6371.3 * np.pi / 180 * bwm
            kdt = scipy.spatial.cKDTree(revents[['xx','yy']].to_numpy()) 
            dists, neighs = kdt.query(revents[['xx','yy']].to_numpy(), nnp+1) # number of nearest neighbors 
            rbwd = list(np.maximum(dists[:,nnp], bwm)) # rbwd = pmax(rbwd, bwm)
        else:
            if len(bwd) != revents.shape[0]:
                raise Exception('length(bwd) != nrow(revents)')
            rbwd = bwd


        # check initial values for the model parameters
        if any([value < 0. for _, value in param0.items()]): # len(list(param0)) != 8 or
            raise Exception("param0 must be a numeric vector of length 8 with positive components")

        param1 = param0.copy()
        thetar = np.empty((no_itr,len(param0.keys()),))
        thetar[:] = np.nan
        asd = np.empty((no_itr,len(param0.keys()),))
        asd[:] = np.nan
        par_names = list(param1.keys())
        loglikfv = np.zeros(revents.shape[0])
        loglikfv[:] = np.nan
        ihess = np.diag(np.full(len(param0.keys()),1.))
        bk = np.zeros(revents.shape[0])
        events = {j[0]: j[1].to_list() for j in revents.iteritems()}
        
        for itr in tqdm(range(0, no_itr)):
            print("declustering:\n")
            for l in range(0, no_itr - itr):
                bkg = decluster(param1, rbwd, events, rpoly, rtperiod, m0,
                                ndiv, model, pb_fix, voronoi)
                events = bkg['revents']
            integ0 = bkg['integ0']
            dbk = bk - np.array(events['bkgd'])
            bk = np.array(events['bkgd'])
            pb = np.array(events['prob'])
            
            # if itr == 0:
            #     save_pickle(param1, 'test/param1')
            #     rdata = dict(revents=events, rpoly=rpoly, tperiod=rtperiod, integ0=integ0)
            #     save_pickle(rdata, 'test/rdata')
            #     save_pickle(ihess, 'test/ihess')
            #     save_pickle(rbwd, 'test/rbwd')
            
            if verbose:
                print("iteration: ", itr)
                print("======================================================")
                print("background seismicity rate:")
                print(pd.DataFrame(bk).describe().transpose())
                print("\nprobability of being a background event:")
                print(pd.DataFrame(pb).describe().transpose())
                print("\nintegral of background seismicity rate: ", integ0)
                print("======================================================")
            
            # if plot_it:
            #     par(mfrow=c(1, 2), mar=c(4, 4, 3, 1))
            #     cols = ifelse(pb < 0.5, "red", "blue")
            #     plot(object$longlat$long, object$longlat$lat, 
            #          cex = 0.05 + 2.5 * revents[, 4]/m0, col=cols,
            #          main=paste("iteration: ", itr), xlab="long", ylab="lat")
            #     polygon(object$region.poly$long, object$region.poly$lat, border=3)
            #     plot(revents[,1], pb, xlab="time",
            #          ylab="probability of being a background event")
            #     rates.inter(param1, object, rbwd, plot.it=plot.it)

            print("estimating:")
            if use_scipy:
                opt = etasfit_scipy(param1, events, rpoly, rtperiod, integ0, m0,
                              ihess, verbose, ndiv, eps, model, impbounds,
                              modified_etas, voronoi)
            else:
                opt = etasfit(param1, events, rpoly, rtperiod, integ0, m0,
                              ihess, verbose, ndiv, eps, model, modified_etas,
                              voronoi)
            
            thetar[itr,:] = opt['estimate']
            loglikfv[itr] = opt['loglik']
            asd[itr,:] = np.sqrt(np.diag(opt['avcov']))
            ihess = np.array(opt['ihessian'])
            gradient = opt['gradient']
            
            if model==1:
                param1 = dict(mu=thetar[itr,0], A=thetar[itr,1], c=thetar[itr,2],
                              alpha=thetar[itr,3], p=thetar[itr,4], D=thetar[itr,5])
            elif model==2:
                param1 = dict(mu=thetar[itr,0], A=thetar[itr,1], c=thetar[itr,2],
                              alpha=thetar[itr,3], p=thetar[itr,4], D=thetar[itr,5])
            elif model==3:
                param1 = dict(mu=thetar[itr,0], A=thetar[itr,1], c=thetar[itr,2],
                              alpha=thetar[itr,3], p=thetar[itr,4], D=thetar[itr,5],
                              q=thetar[itr,6])
                param1 = dict(mu=thetar[itr,0], A=thetar[itr,1], c=thetar[itr,2],
                              alpha=thetar[itr,3], p=thetar[itr,4], D=thetar[itr,5],
                              q=thetar[itr,6])
            else: # model==5 6 and 7
                param1 = dict(mu=thetar[itr,0], A=thetar[itr,1], c=thetar[itr,2],
                              alpha=thetar[itr,3], p=thetar[itr,4], D=thetar[itr,5],
                              q=thetar[itr,6], gamma=thetar[itr,7])

            if verbose:
                print("======================================================")
                print("MLE:\n")
                print(param1)
                print("======================================================")
                # other formulation of etas have the productivity k0 instead of A
                # k0 = 1/np.pi*A*(p-1)*c**(p-1)*(q-1)*D**(q-1)
            
            if itr > 0:
                dtht = ((thetar[itr, :] - thetar[itr - 1, ])/thetar[itr - 1, ]).max()
                dlrv = np.abs(loglikfv[itr] / loglikfv[itr - 1] - 1)
                if np.all(bk == 0.):
                    dbkv = 0.
                else:
                    dbkv = (np.abs(dbk / bk)).max()
                # print([dtht, dlrv, dbkv])
                if all([jj < rel_tol for jj in [dtht, dlrv, dbkv]]):
                    break
                else:
                    if itr == no_itr:
                        print("Reached maximum number of iterations\n")
        
        if verbose:
            print("Execution time:", time.time() - ptm,'s')
        
        # self.catalog.revents = pd.DataFrame(events)
        self.catalog = catalog
        self.b, self.bstd = b, bstd
        self.param1 = param1
        self.revents = events
        self.param = param1
        self.bk = bk
        self.pb = pb
        self.opt = opt
        self.bwd = rbwd
        self.thetar = thetar
        self.loglikfv = loglikfv
        self.asd = asd
        self.integ0 = integ0
        self.gradient = gradient
        self.ndiv = ndiv
        self.itr = itr
        self.exectime = time.time()-ptm
    
    
    def print(self):
    
        # flag = np.array(self.revents['flag'])
        # mm = np.array(self.revents['mm'])[flag == 1]
        # bt = 1. / mm.mean()
        # asd_bt = bt**2 / mm.shape[0]
        try:
            b, bstd = self.b, self.bstd
        except:
            cat = deepcopy(self.catalog)
            b, bstd = self._get_b_value(cat)
        
        param = self.param
        std = {key: self.asd[self.itr,i]  for i, key in enumerate(param.keys())}
        
        param1 = pd.DataFrame(dict(b=[b, bstd]), index=['Estimate', 'StdErr'])
        param2 = pd.DataFrame(self.param, index=['Estimate'])
        param3 = pd.DataFrame(std, index=['StdErr'])
        param23 = pd.concat([param2, param3], axis=0)
        params = pd.concat([param1, param23], axis=1)
        
        string = "ETAS model: fitted using iterative stochastic declustering method\n" + \
                 "converged after " + str(self.itr+1) + " iterations: elapsed execution time " + \
                 str(round(self.exectime/60, 2)) + " minutes\n\n" + \
                 "ML estimates of model parameters:\n" + str(params) + \
                 "\nlog-likelihood: " + str(self.opt['loglik']) + \
                 "\tAIC: " + str(self.opt['aic']) + "\n"
                 # "\nDeclustering probabilities:\n" + str(pd.DataFrame(self.pb).describe().transpose()) + \
        return dict(string=string, params=params, loglik=self.opt['loglik'], aic=self.opt['aic'])


    @classmethod
    def _get_b_value(cls, catalog):
        # b-value of the target events
        flag = np.array(catalog.revents['flag'])
        catalog.purge_catalogue(flag == 1)
        aki_ml = AkiMaxLikelihood()
        b, bstd = aki_ml.calculate(catalog)
        return b, bstd
    


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
            
            
            
#%%        


if __name__ == "__main__":
    
    ###########################################
    # japan
    ###########################################
    
    # df = pd.read_csv('japan_quakes.csv')
    # df['eventID'] = df['Unnamed: 0'].apply(lambda x: str(x))
    # df['date'] = pd.to_datetime(df['date'])
    # df['year'] = df['date'].apply(lambda x: x.year)
    # df['month'] = df['date'].apply(lambda x: x.month)
    # df['day'] = df['date'].apply(lambda x: x.day)
    # df['time'] = pd.to_datetime(df['time'])
    # df['hour'] = df['time'].apply(lambda x: x.hour)
    # df['minute'] = df['time'].apply(lambda x: x.minute)
    # df['second'] = df['time'].apply(lambda x: x.second)
    # df['longitude'] = df['long']
    # df['latitude'] = df['lat']
    # df['magnitude'] = df['mag']
    # df['depth'] = -df['depth']
    
    # # df = df.append(df.iloc[1])
    # # df['year'].iat[-1] += 1
    
    # keys = ['eventID','year', 'month', 'day', 'hour', 'minute', 'second',
    #         'longitude', 'latitude', 'depth', 'magnitude']
    # data = dict()
    # for key in keys:
    #     data[key] = df[key].to_numpy()
    
    # jpoly = dict(long = [134.0, 137.9, 143.1, 144.9, 147.8, 137.8, 137.4, 135.1, 130.6],
    #              lat = [31.9, 33.0, 33.2, 35.2, 41.3, 44.2, 40.2, 38.0, 35.4])

    # catalog = CatalogueEtas(data,
    #                         study_start= datetime.datetime(1953, 5, 26),
    #                         study_end = datetime.datetime(1990, 1, 8),
    #                         region_poly = jpoly, mag_threshold = 4.5)

    # param0 = dict(mu=0.592844590, A=0.204288231, c=0.022692883,
    #               alpha=1.495169224, p=1.109752319, D=0.001175925,
    #               q=1.860044210, gamma=1.041549634)

    # test = Etas(catalog, param0=param0)








    ###########################################
    # iran
    ###########################################

    df = pd.read_csv('test/iran_quakes.csv')
    df['eventID'] = df['Unnamed: 0'].apply(lambda x: str(x))
    df['date'] = pd.to_datetime(df['date'])
    df['year'] = df['date'].apply(lambda x: x.year)
    df['month'] = df['date'].apply(lambda x: x.month)
    df['day'] = df['date'].apply(lambda x: x.day)
    df['time'] = pd.to_datetime(df['time'])
    df['hour'] = df['time'].apply(lambda x: x.hour)
    df['minute'] = df['time'].apply(lambda x: x.minute)
    df['second'] = df['time'].apply(lambda x: x.second)
    df['longitude'] = df['long']
    df['latitude'] = df['lat']
    df['magnitude'] = df['mag']
    
    # df = df.append(df.iloc[1])
    # df['year'].iat[-1] += 1
    
    keys = ['eventID','year', 'month', 'day', 'hour', 'minute', 'second',
            'longitude', 'latitude', 'magnitude']
    data = dict()
    for key in keys:
        data[key] = df[key].to_numpy()
    
    gregion = dict(long = [52., 59., 58., 45., 43.],
                    lat = [26., 25., 29., 38., 35.])
    
    catalog = CatalogueEtas(data,
                            study_start= datetime.datetime(1991, 1, 1),
                            study_end = datetime.datetime(2011, 1, 1),
                            region_poly = gregion, mag_threshold = 4.5)

    param01 = dict(mu=0.5, A=0.2, c=0.05, alpha=2.7, p=1.2,
                    D=0.02, q=2.3, gamma=0.03)

    # etas_iran = Etas(catalog, param0=param01)
    # save_pickle(etas_iran, 'test/etas_iran')

    




    etas_iran = load_pickle('test/etas_iran')
    







    ###########################################
    # italy
    ###########################################

    # df = pd.read_csv('italy_quakes.csv')
    # df['eventID'] = df['Unnamed: 0'].apply(lambda x: str(x))
    # df['date'] = pd.to_datetime(df['date'])
    # df['year'] = df['date'].apply(lambda x: x.year)
    # df['month'] = df['date'].apply(lambda x: x.month)
    # df['day'] = df['date'].apply(lambda x: x.day)
    # df['time'] = pd.to_datetime(df['time'])
    # df['hour'] = df['time'].apply(lambda x: x.hour)
    # df['minute'] = df['time'].apply(lambda x: x.minute)
    # df['second'] = df['time'].apply(lambda x: x.second)
    # df['longitude'] = df['long']
    # df['latitude'] = df['lat']
    # df['magnitude'] = df['mag']
    
    # # df = df.append(df.iloc[1])
    # # df['year'].iat[-1] += 1
    
    # keys = ['eventID','year', 'month', 'day', 'hour', 'minute', 'second',
    #         'longitude', 'latitude', 'magnitude']
    # data = dict()
    # for key in keys:
    #     data[key] = df[key].to_numpy()
    
    # gregion = dict(long = [52., 59., 58., 45., 43.],
    #                lat = [26., 25., 29., 38., 35.])
    
    # catalog = CatalogueEtas(data, dist_unit='km')

    # param01 = dict(mu=0.5, A=0.2, c=0.05, alpha=2.7, p=1.2,
    #                D=0.02, q=2.3, gamma=0.03)

    # test = Etas(catalog, param0=param01)

