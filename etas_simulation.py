# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import os
import sys
sys.path.append('C:\\Users\\Salvatore\\Dropbox\\SalvIac')
import datetime

import numpy as np
import pandas as pd
from scipy.stats import uniform, poisson
from shapely.geometry import Polygon, Point
from tqdm import tqdm

from openquake.hazardlib.scalerel.wc1994 import WC1994
from openquake.hmtk.seismicity.occurrence.aki_maximum_likelihood import AkiMaxLikelihood

from pyrisk.etas.etas8p.catalog import CatalogueEtas
from pyrisk.etas.utils import project
from pyrisk.etas.simulation_functions import (inv_cdf_magnitude_trunc,
                                              inv_cdf_time,
                                              inv_cdf_time_trunc,
                                              inv_cdf_space5,
                                              inv_cdf_space5_trunc,
                                              inv_cdf_space3,
                                              compl_vs_time_p16,
                                              compl_vs_time_hs06,
                                              compl_vs_time_general)
from pyrisk.utils.gardner_knopoff_window import GardnerKnopoffWindowOrig
from myutils.run_multiprocess import run_multiprocess

# only_first_generation
# fully_epidemic
# with or without background



class EtasSimulation():
    
    # simulation options
    allowed_keys = ['num_realization', 'multiprocessing', 'cores',
                    'only_first_generation', 'background', 'save_mode']
    
    
    def __init__(self, fit, catalog_to_start, model={}, mag_max=7.5,
                 sim_start=None, sim_end=None,
                 simulaton_region=None, buffer_region=None,
                 folder=None, simul_options={}):
        
        # checks
        if ((sim_start is None) and (sim_end is not None)) or \
           ((sim_start is not None) and (sim_end is None)):
            raise Exception('check sim_start and sim_end!')
        
        
        # etas model
        self.fit = fit
        self.fit_catalog = fit.catalog
        self.bwd = fit.bwd
        self.param = fit.param
        self.mag_threshold = fit.catalog.mag_threshold

        # etas model
        if model is None:
            model = {}
        model.setdefault('time', 'trunc') # 1 truncated 2 untruncated
        model.setdefault('space', 5) # numbers as Zhuang et al. 2011
        model.setdefault('magnitude', 'trunc') # 1 truncated 2 untruncated #TODO
        self.model = model
        # for magnitude model
        self.mag_max = mag_max
        self.b = self._get_b_value(fit)

        
        # simulation time window
        self.sim_start = sim_start
        self.sim_end = sim_end
        self.t = CatalogueEtas.date2day(sim_start, sim_start)
        self.dt = CatalogueEtas.date2day(sim_end, sim_start)
        
        
        # options (set default and override)
        if simul_options is None:
            simul_options = {}
        simul_options.setdefault('num_realization', 10000)
        simul_options.setdefault('multiprocessing', True)
        simul_options.setdefault('cores', 3)
        simul_options.setdefault('only_first_generation', False)
        simul_options.setdefault('background', True)
        simul_options.setdefault('save_mode', True)
        
        self.simul_options = simul_options     
        self.__dict__.update((k, v) for k, v in simul_options.items() 
                             if k in self.allowed_keys)
        
        if folder is None:
            self.output_folder = os.path.join(os.getcwd())
        else:
            self.output_folder = os.path.join(os.getcwd(), folder)
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        
        # simulation and buffer region
        self.simulation_region = simulaton_region
        self.simulation_region_xy = self._lonlat2xy(simulaton_region['lon'],
                                                    simulaton_region['lat'])
        if buffer_region is None:
            self.buffer_region = simulaton_region
        else:
            self.buffer_region = buffer_region
        self.buffer_region_xy = self._lonlat2xy(self.buffer_region['lon'],
                                                self.buffer_region['lat'])
        
        # preprocess input catalog
        catalog_to_start['tt'] = CatalogueEtas.date2day(catalog_to_start['datetime'], self.sim_start)
        event_x, event_y = project(catalog_to_start['longitude'], catalog_to_start['latitude'], self.fit)        
        catalog_to_start['xx'] = event_x
        catalog_to_start['yy'] = event_y
        catalog_to_start['mm'] = catalog_to_start['magnitude'].to_numpy()-self.fit.catalog.mag_threshold
        # convert fault geomentry (where available)
        geom = [None]*catalog_to_start.shape[0]
        for i, (_, row) in enumerate(catalog_to_start.iterrows()):
            if row["geometry"] is not None:
                proj = CatalogueEtas.longlat2xy(row["geometry"].lons.flatten(),
                                                row["geometry"].lats.flatten(),
                                                region_poly=self.fit.catalog.region_poly,
                                                dist_unit=self.fit.catalog.dist_unit)
                geom[i] = {'x': proj['x'],
                           'y': proj['y'],
                           'depth': row["geometry"].depths}         
        catalog_to_start["geom"] = geom
        self.catalog_to_start = catalog_to_start
        
    
    
    def _define_inputs(self):
        # input for simulations
        self.__dict__.update((k, v) for k, v in self.simul_options.items() 
                             if k in self.allowed_keys)
        self.inputs = list()
        for l in range(0, self.num_realization):
            dic = {'gen0': self.catalog_to_start,
                   'fit': self.fit,
                   'params': self.param,
                   'model': self.model,
                   't': self.t,
                   'tdt': self.t+self.dt,
                   'b': self.b,
                   'mag_max': self.mag_max,
                   'simulation_region_xy': self.simulation_region_xy,
                   'buffer_region_xy': self.buffer_region_xy,
                   'include_background': self.simul_options['background'],
                   'only_first_generation': self.simul_options['only_first_generation'],
                   'output_folder': self.output_folder,
                   'iter': l,
                   "save_mode": self.simul_options['save_mode']}
            self.inputs.append(dic)
        
    

    def simulate(self):
        self._define_inputs()
        if self.multiprocessing:
            out = run_multiprocess(self.simulate_multi, self.inputs, self.cores)
        else:
            out = list()
            for inpu in tqdm(self.inputs):
                out.append( self.simulate_multi(inpu) )
        return out
    

        
    def _lonlat2xy(self, lons, lats):
        proj = CatalogueEtas.longlat2xy(long=np.array(lons),
                                        lat=np.array(lats),
                                        region_poly=self.fit.catalog.region_poly,
                                        dist_unit=self.fit.catalog.dist_unit)
        coords = [(x, y) for x,y in zip(proj['x'], proj['y'])]
        region_win = Polygon(coords)       
        return region_win
        
        
    def _get_b_value(self, fit):
        # b-value of the target events
        flag = np.array(fit.catalog.revents['flag'])
        fit.catalog.purge_catalogue(flag == 1)
        aki_ml = AkiMaxLikelihood()
        b, _ = aki_ml.calculate(fit.catalog)
        return b     
    
        
    def __str__(self):
        string = "<EtasSimulation"+" \n"\
            "sim_start: "+str(self.sim_start)+", sim_end: "+str(self.sim_end)+"\n"\
            "model: "+str(self.model)+"\n"+\
            "max_mag: "+str(self.mag_max)+", b-value: "+str(self.b)+" \n"+\
            "ETAS params: \n"+str(self.param)+"\n"+\
            "simulation options: \n"+str(self.simul_options)+">"
        return string
    
    
    def __repr__(self):
        return self.__str__()        


    
    
    # simulate
    @staticmethod
    def simulate_multi(inputs):
    
        gen0 = inputs['gen0']
        fit = inputs['fit']
        params = inputs['params']
        model = inputs['model']
        t = inputs['t']
        tdt = inputs['tdt']
        b = inputs['b']
        mag_max = inputs['mag_max']
        mag_threshold = fit.catalog.mag_threshold
        simulation_region_xy = inputs['simulation_region_xy']
        buffer_region_xy = inputs['buffer_region_xy']
        include_background = inputs['include_background']
        only_first_generation = inputs['only_first_generation']
        path = inputs['output_folder']
        it = inputs['iter']
        save_mode = inputs['save_mode']
        
        # for reproducibility
        import numpy as np
        np.random.seed(seed=it)
        
        if include_background:
            bkgr = EtasSimulation.get_background(fit, t, tdt, buffer_region_xy)
            gen0 = pd.concat([gen0, bkgr], ignore_index=True)
            
        gg = [gen0] # list with each generation, regardless the parent event
        
        timedep_mc = None
        if model["magnitude"]["type"] == "trunc_complete_advanced":
            timedep_mc = EtasSimulation.calc_timedep_mc(gg[-1], gg[-1], mag_threshold, 5)
            
        gl = EtasSimulation.get_following_generation(gg[-1], params, model, t, tdt,
                                      simulation_region_xy, buffer_region_xy,
                                      b, mag_max, mag_threshold, timedep_mc)
        gg.append(pd.concat(gl, ignore_index=True))
        
        if not only_first_generation:
            l = 1
            while len(gl) != 0: # all([gli.empty for gli in gl]):
                #print(l)
                if model["magnitude"]["type"] == "trunc_complete_advanced":
                    timedep_mc = EtasSimulation.calc_timedep_mc(gg[-1],
                                                 pd.concat(gg, ignore_index=True),
                                                 mag_threshold, 0)  
                    
                gl = EtasSimulation.get_following_generation(gg[-1], params, model, t, tdt,
                                              simulation_region_xy, buffer_region_xy,
                                              b, mag_max, mag_threshold, timedep_mc)
                if len(gl) != 0:
                    gg.append(pd.concat(gl, ignore_index=True))
                    l += 1
    
        stoch_catalog = pd.concat(gg, ignore_index=True)
        stoch_catalog.drop_duplicates(subset=['tt', 'xx', 'yy', 'mm'], inplace=True)
        stoch_catalog.sort_values(by=['tt'], inplace=True)
        
        # convert tt (days) in datetime
        stoch_catalog["datetime"] = stoch_catalog["datetime"].iloc[0] - \
                                    datetime.timedelta(days=stoch_catalog["tt"].iloc[0]) + \
                                    np.array([datetime.timedelta(days=days) for days in stoch_catalog["tt"]])
        
        # filter
        stoch_catalog = EtasSimulation.filter_events(stoch_catalog, t, tdt, simulation_region_xy, mag_threshold)
        
        # convert x y in lon lat
        proj = CatalogueEtas.xy2longlat(stoch_catalog['xx'],
                                        stoch_catalog['yy'],
                                        region_poly=fit.catalog.region_poly,
                                        dist_unit=fit.catalog.dist_unit)
        stoch_catalog['longitude'] = proj['long']
        stoch_catalog['latitude'] = proj['lat']    
        stoch_catalog.reset_index(inplace=True, drop=True)
        
        if save_mode: # save catalog
            stoch_catalog.to_csv(path+'\stoch_catalog_'+str(it).zfill(5)+'.csv')
            return it
        else:
            return stoch_catalog
    
    
    
    @staticmethod
    def get_background(fit, t_start, t_end, region_win_xy):
        '''
        Generate the background catalog with the estimated background intensity
        μ(x,y) (equation 11 in Zhuang et al. 2011)
        For each event in the background catalog, generate a random variable Ui 
        uniformly distributed in [0,1], accept it if Ui < ν φi dt/(t−t0), 
        where ν is as defined in eq (12) Zhuang et al. 2011
        t0 is the starting time of the catalog and t−t0 is the length
        of period of data fitted to the model.
        Randomly assign each selected event a new time uniformly distributed
        in [t,t+dt], and relocate each selected even by adding a 2D Gaussian
        random variable with a density Zdi, where Z is the kernel function
        used in estimating the background seismicity and di is the bandwidth
        corresponding to the selected event
        '''
        
        # get info from catalog and fit
        cat = fit.catalog.revents[fit.catalog.revents['flag'] == 1]
        pb = fit.pb[fit.catalog.revents['flag'] == 1]
        bwd = np.array(fit.bwd)[fit.catalog.revents['flag'] == 1]
        tt0 = fit.catalog.rtperiod['study_end'] - fit.catalog.rtperiod['study_start']
        
        # pick random events
        u = uniform.rvs(size=cat.shape[0])
        higher = u < fit.param['mu'] * pb * (t_end-t_start)/(tt0)
        
        # add 2D gaussian (https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf)
        x_old = cat['xx'].iloc[higher].to_numpy()
        y_old = cat['yy'].iloc[higher].to_numpy()
        bwd = bwd[higher]
        theta = uniform.rvs(loc=0., scale=2*np.pi, size=bwd.shape[0])
        u = uniform.rvs(size=bwd.shape[0])
        r = np.sqrt( -2 * np.log(1-u) ) # inverse cdf of the relative locations (f)
        event_x = x_old + r*bwd*np.cos(theta)
        event_y = y_old + r*bwd*np.sin(theta)
        
        # filter for region win
        filt = np.array([Point(xxx, yyy).within(region_win_xy) for xxx, yyy in zip(event_x, event_y)])
        
        if np.sum(filt) == 0:
            backg = pd.DataFrame([], columns = ['tt','xx','yy','mm','magnitude','geometry','geom'])
        else:
            # random times
            times = uniform.rvs(loc=t_start, scale=(t_end-t_start), size=np.sum(filt))
            
            backg = pd.DataFrame(zip(times,
                                     event_x[filt],
                                     event_y[filt],
                                     cat['mm'].to_numpy()[higher][filt],
                                     cat['mm'].to_numpy()[higher][filt] + fit.catalog.mag_threshold,
                                     [None]*np.sum(filt),[None]*np.sum(filt)),
                             columns = ['tt','xx','yy','mm','magnitude','geometry','geom'])
            backg.sort_values(by=['tt'], inplace=True)
            backg.reset_index(inplace=True, drop=True)
        
            # heatmap, xedges, yedges = np.histogram2d(testx, testy, bins=50)
            # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            # plt.clf()
            # plt.imshow(heatmap.T, extent=extent, origin='lower')
            # plt.show()
        
        return backg
        
        
       
    @staticmethod
    def get_following_generation(gen, params, model, t_start, t_end,
                                 simulation_region_xy, buffer_region_xy,
                                 b, mag_max, mag_threshold, timedep_mc):
        ''' For each event i , namely (ti, xi, yi, mi), in the catalog G(l), 
        simulate its N(i) offspring, namely, Oi(l) where N(i) is a Poisson random
        variable with a mean of κ(mi), and tk(i) (xk(i), yk(i)) and mk(i) are
        generated from the probability densities g, f and s respectively. 
        Let Oi(l) ← {(tk , xk , yk , mk)} where tk in [t, t + dt].
        
        # for km (productivity): the expectation of the number of children
        # (unit: events) spawned by an event of magnitude m
        A = fit.param['A']
        alpha = fit.param['alpha']
        
        # for the pdf of the length of the time interval between a child and its
        # parent (unit: day−1)
        c = fit.param['c']
        p = fit.param['p']
        
        # for the pdf of the relative locations between the parent and children
        # (unit: deg−2)
        D = fit.param['D']
        q = fit.param['q']
        gamma = fit.param['gamma']
        '''
        
        if ((t_start is None) and (t_end is not None)) or \
            ((t_start is not None) and (t_end is None)):
            raise Exception('check t_start and t_end!')
    
        # account for full fault geometry
        mag_parent = list()
        mag_to_use = list()
        time_parent = list()
        x_par = list()
        y_par = list()
        points = list()
        timedep_mc_par = list()
        wc = WC1994()
        for par, geom in gen['geom'].items():
            if geom is None:
                mag_parent.append(gen['magnitude'][par])
                mag_to_use.append(gen['magnitude'][par])
                time_parent.append(gen['tt'][par])
                x_par.append(gen['xx'][par])
                y_par.append(gen['yy'][par])
                points.append(1)
                if timedep_mc is not None:
                    timedep_mc_par.append(timedep_mc[par])
                else:
                    timedep_mc_par.append(None)
            else:
                for i in range(0, geom['x'].shape[0]):
                    mag_parent.append(gen['magnitude'][par])
                    time_parent.append(gen['tt'][par])
                    x_par.append(geom['x'][i])
                    y_par.append(geom['y'][i])
                    points.append(geom['x'].shape[0])
                    if timedep_mc is not None:
                        timedep_mc_par.append(timedep_mc[par])
                    else:
                        timedep_mc_par.append(None)
                    # mag_to_use.append(mag_parent[-1])
                    mag_to_use.append( np.max([ mag_threshold,
                                       wc.get_median_mag(wc.get_median_area(mag_parent[-1], None)/points[-1], None) ]))
                    
        # # productivity (corrected by the number of points UCERF3)
        # km = (params["A"] / np.array(points)) * np.exp(params["alpha"] * \
        #                                                (np.array(mag_parent)-mag_threshold))
        if timedep_mc is None:
            km = (params["A"] / np.array(points)) * np.exp(params["alpha"] * \
                                              (np.array(mag_parent)-mag_threshold))
        else:
            km = (params["A"] / np.array(points)) * np.exp(params["alpha"] * \
                                              (np.array(mag_parent)-np.array(timedep_mc_par)))
    
        # km[np.array(mag_parent) <= 5.49] = 0. 
        # this assumes that the mainshock occurs at zero (time_parent)
        # mag_main = model["magnitude"]["mag_main"]
        # ind = np.array(time_parent) > 0.
        # mc = np.clip(mag_main-4.5-0.75*np.log10(np.array(time_parent)[ind]),
        #              mag_threshold, np.array(mag_parent)[ind])
        # km[ind] = (params["A"] / np.array(points)[ind]) * \
        #                 np.exp(params["alpha"] * (np.array(mag_parent)[ind]-mc))
        
        
        
        # random number of offspring events
        ni = poisson.rvs(km)
        if isinstance(ni, int):
            ni = np.array([ni])
    
    
        ol = list() # list of offsprings for each event of the generation l
        for par in range(0, ni.shape[0]):
            if ni[par] > 0:
                o = EtasSimulation.generate_offsprings(params, model, ni[par], time_parent[par],
                                        x_par[par], y_par[par], mag_to_use[par],
                                        mag_parent[par], b, mag_max, mag_threshold)
                                        # timedep_mc_par[par])
                o = EtasSimulation.filter_events(o, t_start, t_end, buffer_region_xy, mag_threshold)
                # o.sort_values(by=['tt'], inplace=True)
                # o.reset_index(inplace=True, drop=True)
                ol.append(o)
    
        return ol
    
    
    
    @staticmethod
    def filter_events(o, t_start, t_end, region_win_xy, mag_threshold):
        if (t_start is not None) and \
            (t_end is not None) and \
            (region_win_xy is not None):
            # filter time-window and region
            o = o[ (o['tt'] >= t_start) & (o['tt'] <= t_end) & 
                    np.array([Point(xxx, yyy).within(region_win_xy) for xxx, yyy in zip(o['xx'], o['yy'])]) ]
        elif region_win_xy is not None:
            o = o[ np.array([Point(xxx, yyy).within(region_win_xy) for xxx, yyy in zip(o['xx'], o['yy'])]) ]
        elif (t_start is not None) and (t_end is not None):
            o = o[ (o['tt'] >= t_start) & (o['tt'] <= t_end) ]
        o = o[ o['magnitude']>=mag_threshold ]
        return o
    
    
    
    @staticmethod
    def generate_offsprings(params, model, num_offspr, time_parent, x_par, y_par,
                            mag_to_use, mag_parent, b, mag_max, mag_threshold):
        
    
        # g = (p - 1)/c * np.power(1. + ttt/c, -p)
        u = uniform.rvs(size=num_offspr)
        # inverse cdf of time interval pdf (g)
        # deltat = params["c"]*((1-u)**(1/(1-params["p"]))-1)
        if model["time"] == "trunc":
            # gk = GardnerKnopoffWindowOrig()
            # ta = gk.calc(m_offspr)[1]*364.75
            ta = 5*365
            deltat = inv_cdf_time_trunc(u, params["c"], params["p"], ta)
        elif model["time"] == "untrunc":
            deltat = inv_cdf_time(u, params["c"], params["p"])
        t_offspr = time_parent + deltat
        
        
        # only account for t_offspr after 0. (to speed things up)
        #TODO this does not work if history is not provided until the start of the simulation period
        deltat = deltat[t_offspr >= 0.]
        t_offspr = t_offspr[t_offspr >= 0.]
        num_offspr = t_offspr.shape[0]
        
        
        if num_offspr > 0:
            # b = (1. / gl['mm'].mean())/np.log(10.)
            mmax = mag_max # np.random.choice([mag_parent, mag_max])
            mmin = mag_threshold
            # m_offspr = mmin + (np.log10(1.-u*(1.-10.**(-b*(mmax-mmin)))))/(-b)
            # if isinstance(model["magnitude"], str): # use functions
            #     u = uniform.rvs(size=num_offspr)
            #     if model["magnitude"] == "trunc":
            #         m_offspr = inv_cdf_magnitude_trunc(u, b, mmin, mmax)
            #     elif model["magnitude"] == "trunc_complete":
            #         mc = np.clip(compl_vs_time_hs06(mag_parent, deltat), mmin, mmax-1.)
            #         m_offspr = inv_cdf_magnitude_trunc(u, b, mc, mmax) 
    
            if isinstance(model["magnitude"], dict): # use functions
                u = uniform.rvs(size=num_offspr)
                if model["magnitude"]["type"] == "trunc":
                    m_offspr = inv_cdf_magnitude_trunc(u, b, mmin, mmax)
                elif model["magnitude"]["type"] == "trunc_complete_advanced":
                    # print("ok")
                    if mag_parent >= 6:
                        mc = np.clip(compl_vs_time_hs06(mag_parent, deltat), mmin, mag_parent)
                    else:
                        mc = mmin
                    m_offspr = inv_cdf_magnitude_trunc(u, b, mc, mmax)
            # if isinstance(model["magnitude"], dict): # use functions
            #     u = uniform.rvs(size=num_offspr)
            #     if model["magnitude"][0] == "trunc_complete":
            #         mag_mainshock = model["magnitude"][1]
            #         coeffs = model["magnitude"][2]
            #         mc = np.clip(compl_vs_time_general(mag_mainshock, t_offspr, *coeffs), mmin, mmax-.1)
            #         m_offspr = inv_cdf_magnitude_trunc(u, b, mc, mmax)    
                 
                elif model["magnitude"] == "untrunc":
                    raise Exception("mag model untrunc not implemented")
                else:
                    raise Exception("unknown mag model")
            
            elif isinstance(model["magnitude"], list): 
                # # this portion of the code does not work with t_offspr < 0 #TODO
                # # if np.any(t_offspr < 0.):
                # #     raise Exception("blabla")
                # u = uniform.rvs(size=num_offspr)
                # # mc = np.clip(compl_vs_time_general(mag_parent, t_offspr, *[0.5, 0.25, 1.]),
                # #               mmin, mmax-0.01)
                # mc = mmin
                # m_offspr = np.zeros(num_offspr)
                
                # m_offspr[t_offspr<=1.] = inv_cdf_magnitude_trunc(u[t_offspr<=1.],
                #                                    model["magnitude"][0],
                #                                    mc, mmax) 
                # m_offspr[t_offspr>1.] = inv_cdf_magnitude_trunc(u[t_offspr>1.], 
                #                                    model["magnitude"][1], 
                #                                    mc, mmax) 
                if mag_parent >= 6:
                    m_offspr = np.zeros(num_offspr)
                    mfdt = model["magnitude"][2]
                    timecut = model["magnitude"][1]
                    m_offspr[t_offspr<=timecut] = np.clip(mag_parent + \
                                                   np.random.choice(mfdt, size=np.sum(t_offspr<=timecut)),
                                                   mmin, mmax)
                    u = uniform.rvs(size=np.sum(t_offspr>timecut))
                    m_offspr[t_offspr>timecut] = inv_cdf_magnitude_trunc(u, 
                                                                model["magnitude"][0], 
                                                                mmin, mmax) 
                else:
                    u = uniform.rvs(size=num_offspr)
                    m_offspr = inv_cdf_magnitude_trunc(u, model["magnitude"][0], 
                                                        mmin, mmax) 
    
        
            
            
            else: # use random sampling from array
                u = uniform.rvs(size=np.sum(t_offspr<=0.))
                m_offspr = inv_cdf_magnitude_trunc(u, b, mmin, mmax)
                if np.any(t_offspr>0.):
                    m_offspr = list(m_offspr)
                    mfdt = model["magnitude"]
                    for dt in t_offspr[t_offspr>0.]:
                        ind = np.where(mfdt["time"] < dt)[0][-1]
                        if mfdt["magnitude"][ind].shape[0] != 0:
                            m_offspr.append( np.random.choice(mfdt["magnitude"][ind], size=1)[0] )
                        else:
                            m_offspr.append( np.random.choice(mfdt["magnitude"][ind+1], size=1)[0] )
                        # shift = test_shift(mfdt["time"], dt, time_parent, 5)
                        # ind = np.logical_and(mfdt["time"] <= dt+shift, mfdt["time"] >= time_parent-shift)
                        # m_offspr.append( np.random.choice(mfdt["magnitude"][ind], size=1)[0] )
                    m_offspr = np.array(m_offspr)
            
            
            u_r = uniform.rvs(size=num_offspr)
            u_theta = uniform.rvs(size=num_offspr) # loc=0., scale=2*np.pi, size=num_offspr)
            dm = mag_to_use-mag_threshold
            if model["space"] == 1:
                raise Exception("not implemented")
            elif model["space"] == 2:
                raise Exception("not implemented")
            elif model["space"] == 3:
                # f = (q - 1) / (D * np.exp(gamma * dm) * np.pi) * np.power(1 + r**2 / (D * np.exp(gamma * dm)), - q)
                # inverse cdf of the relative locations (f)
                # r = np.sqrt( params["D"]*((1-u)**(-1/(params["q"]-1))-1)/np.exp(-params["gamma"]*dm) )
                xsim, ysim = inv_cdf_space3(u_r, u_theta, dm, params["q"], params["D"])
            elif model["space"] == 4:
                raise Exception("not implemented")
            elif model["space"] == 5:
                # f = (q - 1) / (D * np.exp(gamma * dm) * np.pi) * np.power(1 + r**2 / (D * np.exp(gamma * dm)), - q)
                # inverse cdf of the relative locations (f)
                # r = np.sqrt( params["D"]*((1-u)**(-1/(params["q"]-1))-1)/np.exp(-params["gamma"]*dm) )
                # xsim, ysim = inv_cdf_space5(u_r, u_theta, dm, params["q"], params["D"], params["gamma"])
                r_trunc = 1.
                xsim, ysim = inv_cdf_space5_trunc(u_r, u_theta, dm, params["q"],
                                                  params["D"], params["gamma"], r_trunc)
            x_offspr = x_par + xsim
            y_offspr = y_par + ysim
        
        else:
            x_offspr = np.array([])
            y_offspr = np.array([])
            m_offspr = np.array([])
        
        # out
        o = pd.DataFrame(zip(t_offspr, x_offspr, y_offspr, m_offspr-mag_threshold,
                             m_offspr, [None]*t_offspr.shape[0], [None]*t_offspr.shape[0]),
                         columns = ['tt','xx','yy','mm','magnitude','geometry','geom'])
        return o
    
    
    
    # def test_shift(mfdt_time, time, time_parent, shift):
    #     print(shift)
    #     if np.sum((mfdt_time <= time+shift) & (mfdt_time >= time_parent-shift)) == 0:
    #         if shift == 0.:
    #             shift = 0.1
    #         shift = test_shift(mfdt_time, time, time_parent, shift*2)
    #     return shift
    
    @staticmethod
    def calc_timedep_mc(gg_next_generation, gg_past_events, mag_threshold, indtt):
        mcs = list()
        for g in gg_next_generation.itertuples(index=False): #.iterrows(): # 
            ind = (gg_past_events["tt"] < g[indtt]) & \
                  (g[indtt]-gg_past_events["tt"] < 1.) & \
                  (gg_past_events["magnitude"] >= 6.)
            # ind = (gg_past_events["tt"] < g["tt"]) & \
            #       (g["tt"]-gg_past_events["tt"] < 1.) & \
            #       (gg_past_events["magnitude"] >= 6.)
            if np.any(ind):
                time = g[indtt] - gg_past_events.loc[ind]["tt"]
                mag_par = gg_past_events.loc[ind]["magnitude"]
                mc = mag_par-4.5-0.75*np.log10(time)
                mc = np.clip(mc, mag_threshold, mag_par)
                mc_max = np.max(mc)
                mcs.append(mc_max)
            else:
                mcs.append(mag_threshold)
        return pd.Series(mcs, index=gg_next_generation.index)
    
    # @staticmethod
    # def calc_timedep_mc(gg_next_generation, gg_past_events, mag_threshold, indtt):
    #     mcs = list()
    #     for _, g in gg_next_generation.iterrows(): 
    #         ind = (gg_past_events["tt"] < g["tt"]) & \
    #               (g["tt"]-gg_past_events["tt"] < 1.) & \
    #               (gg_past_events["magnitude"] >= 6.)
    #         if np.any(ind):
    #             time = g["tt"] - gg_past_events.loc[ind]["tt"]
    #             mag_par = gg_past_events.loc[ind]["magnitude"]
    #             mc = mag_par-4.5-0.75*np.log10(time)
    #             mc = np.clip(mc, mag_threshold, mag_par)
    #             mc_max = np.max(mc)
    #             mcs.append(mc_max)
    #         else:
    #             mcs.append(mag_threshold)
    #     return pd.Series(mcs, index=gg_next_generation.index)


