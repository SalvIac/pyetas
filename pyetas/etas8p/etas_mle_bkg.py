# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""

import time
from copy import deepcopy

import numpy as np
import pandas as pd
# from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from shapely.geometry import Point

from pyrisk.etas.etas8p.dist import dist
from pyrisk.etas.etas8p.etasfit import etasfit
from pyrisk.etas.etas8p.etasfit_scipy import etasfit_scipy # q=1.5 alpha=beta
from pyrisk.etas.etas8p.voronoi import get_voronoi
from pyrisk.etas.etas8p.etas import Etas
from myutils.utils_pickle import load_pickle, save_pickle


class EtasMleBkg(Etas):
        
    def __init__(self, catalog, param0=None, bkg=None, verbose=True,
                 pb_fix=False, ndiv=1000, rel_tol=1e-03, eps=1e-06,
                 use_scipy=False, model=6, fault=True, impbounds=None,
                 voronoi_mode=False, modified_etas=False):
        
        if bkg is None:
            raise Exception("bkg must be specified to use EtasMleBkg!")
        if bkg['rated'].shape[0] != bkg['lon'].shape[0] or bkg['rated'].shape[1] != bkg['lat'].shape[0]:
            raise Exception("bkg[rated] must be specified with lons on rows and lats on columns!")
        
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

        # check initial values for the model parameters
        if any([value < 0. for _, value in param0.items()]): # len(list(param0)) != 8 or
            raise Exception("param0 must be a numeric vector of length 8 with positive components")

        param1 = param0.copy()
        thetar = np.empty((1,len(param0.keys()),))
        thetar[:] = np.nan
        asd = np.empty((1,len(param0.keys()),))
        asd[:] = np.nan
        loglikfv = np.zeros(revents.shape[0])
        loglikfv[:] = np.nan
        ihess = np.diag(np.full(len(param0.keys()),1.))
        events = {j[0]: j[1].to_list() for j in revents.iteritems()}
    
    
        
        # background calculations
        yy, xx = np.meshgrid(bkg['lat'], bkg['lon'])
        xy = catalog.longlat2xy(long=xx,
                                lat=yy,
                                region_poly=catalog.region_poly,
                                dist_unit=catalog.dist_unit)
        x = xy['x'][:,0]
        y = xy['y'][0,:]
        bkg = {'lon': bkg['lon'],
               'lat': bkg['lat'],
                'x': x,
                'y': y,
                'rated': bkg['rated']}

        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # aa = plt.scatter(xx.flatten(), yy.flatten(),
        #                 c=bkg['rated'].flatten(), cmap='BuPu')
        # plt.show()

        # fig = plt.figure()
        # aa = plt.scatter(xy['x'].flatten(), xy['y'].flatten(),
        #                 c=bkg['rated'].flatten(), cmap='BuPu')
        # plt.show()
        # stop

        
        integ0 = EtasMleBkg.calc_integ_bkg(bkg, win, rtperiod)
        bk = EtasMleBkg.calc_bkg(bkg, catalog)
        events['bkgd'] = bk
    
        if verbose:
            print("======================================================")
            print("background seismicity rate:")
            print(pd.DataFrame(bk).describe().transpose())
            print("\nintegral of background seismicity rate: ", integ0)
            print("======================================================")

        print("estimating:")
        for i in range(0,2): # to be sure it doesn't get stuck
            if use_scipy:
                opt = etasfit_scipy(param1, events, rpoly, rtperiod, integ0, m0,
                              ihess, verbose, ndiv, eps, model, impbounds, 
                              modified_etas, voronoi)
            else:
                opt = etasfit(param1, events, rpoly, rtperiod, integ0, m0,
                              ihess, verbose, ndiv, eps, model, modified_etas,
                              voronoi)
            
            thetar[0,:] = opt['estimate']
            loglikfv[0] = opt['loglik']
            asd[0,:] = np.sqrt(np.diag(opt['avcov']))
            ihess = np.array(opt['ihessian'])
            gradient = opt['gradient']
            
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
            else: # model==5 6 and 7
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
        self.b, self.bstd = b, bstd
        self.param1 = param1
        self.revents = events
        self.param = param1
        self.bk = bk
        self.opt = opt
        self.bwd = None
        self.thetar = thetar
        self.loglikfv = loglikfv
        self.asd = asd
        self.integ0 = integ0
        self.gradient = gradient
        self.ndiv = ndiv
        self.itr = None
        self.exectime = time.time()-ptm
    

    @classmethod
    def calc_integ_bkg(cls, bkg, region_win, tperiod):
        '''
        x and y coordinates (relative lat lon) for bkg['x'] bkg['y']
        region_win in relative lat lon
        time period in days
        rate in bkg['rated'] is a daily rate (per unit area)
        '''
        yy, xx = np.meshgrid(bkg['y'], bkg['x'])
        _sum = 0.
        for i in range(0,bkg['rated'].shape[0]):
            for j in range(0,bkg['rated'].shape[1]):
                if Point(xx[i,j], yy[i,j]).within(region_win):
                    _sum += (bkg['rated'][i,j])
        # the rate is per unit area
        dx = (np.diff([np.min(bkg['x']), np.max(bkg['x'])]) / (bkg['x'].shape[0]-1))[0]
        dy = (np.diff([np.min(bkg['y']), np.max(bkg['y'])]) / (bkg['y'].shape[0]-1))[0]
        return (_sum * dx * dy * (tperiod['study_end']-tperiod['study_start']))
    
    
    @classmethod
    def calc_bkg(cls, bkg, catalog):
        # interpolant = interp2d(r['x'], r['y'], r['bkgd'].transpose())
        intrp = RectBivariateSpline(bkg['x'], bkg['y'], bkg['rated'])
        bk = intrp.ev(catalog.revents['xx'], catalog.revents['yy'])

        # # test time
        
        # tt = time.time()
        # test1 = intrp.ev(fit.catalog.data['longitude'], fit.catalog.data['latitude'])
        # print(time.time()-tt)
    
        # tt = time.time()
        # test3 = list()
        # for lons, lats in zip(fit.catalog.data['longitude'], fit.catalog.data['latitude']):
        #     # if min(lon) <= lons <= max(lon) and min(lat) <= lats <= max(lat):
        #     test3.append(float(intrp.ev(lons, lats)))
        # print(time.time()-tt)
    
        # tt = time.time()
        # test4 = list()
        # for lons, lats in zip(fit.catalog.data['longitude'], fit.catalog.data['latitude']):
        #     # if min(lon) <= lons <= max(lon) and min(lat) <= lats <= max(lat):
        #     test4.append(float(interpolant(lons, lats)))
        # print(time.time()-tt)
        
        # fig = plt.figure()
        # plt.plot((np.array(test3) - np.array(test4))/np.array(test4))
        # plt.show()
        # np.argmin((np.array(test3) - np.array(test4))/np.array(test4))
        
        return bk

