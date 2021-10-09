# -*- coding: utf-8 -*-
"""
@author: Salvatore
https://github.com/jalilian/ETAS/tree/master/R
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ioff()
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from pyetas.etas8p.etas import Etas
from pyetas.etas8p.catalog import CatalogueEtas
from pyetas.etas8p.dist import dist, dist2
from pyetas.etas8p.poly import dGauss
from myutils.utils_pickle import load_pickle, save_pickle




def probs(fit):
    # spatstat::verifyclass(fit, "etas")
    xx = fit.revents['xx']
    yy = fit.revents['yy']
    # probability of each event being a triggered event
    pb = 1. - np.array(fit.revents['prob'])
    target = np.array(fit.revents['flag']) == 1
    return pd.DataFrame(zip(xx,yy,pb,target),
                        columns=['long', 'lat', 'prob', 'target'])


def rates(fit, lat_range=None, long_range=None, dimyx=None, plot_it=True):
    if not isinstance(fit, Etas):
        raise Exception('variable fit is not an Etas instance')
    out = rates_inter(fit.param, fit, fit.bwd, lat_range=lat_range,
                      long_range=long_range, dimyx=dimyx)
    return out


def rates_inter(theta, fit, bwd, lat_range=None, long_range=None, dimyx=None):
    catalog = fit.catalog
    if lat_range is None:
        lat_range = [min(catalog.region_poly['lat']),
                     max(catalog.region_poly['lat'])]
    if long_range is None:
        long_range = [min(catalog.region_poly['lon']),
                      max(catalog.region_poly['lon'])]
    lat_range.extend([lat_range[0], lat_range[1]])
    long_range.extend([long_range[1], long_range[0]])
    lat_range = np.array(lat_range)
    long_range = np.array(long_range)
    xy_bnd = CatalogueEtas.longlat2xy(long=long_range,
                                      lat=lat_range,
                                      region_poly=catalog.region_poly,
                                      dist_unit=catalog.dist_unit)
    if dimyx is None:
        rangex = [np.min(xy_bnd['x']), np.max(xy_bnd['x'])]
        rangey = [np.min(xy_bnd['y']), np.max(xy_bnd['y'])]
        rv = np.diff(rangex)[0] / np.diff(rangey)[0]
        if (rv > 1):
            dimyx = [128., round(128.*rv)]
        else:
            dimyx = [round(128.*(1 / rv)), 128.]
    
    gx = np.linspace(min(xy_bnd['x']), max(xy_bnd['x']), int(dimyx[1]))
    gy = np.linspace(min(xy_bnd['y']), max(xy_bnd['y']), int(dimyx[0]))
    out = cxxrates(theta, fit.revents, bwd, catalog.rtperiod, gx, gy)
    
    out = dict(x=np.linspace(long_range[0], long_range[1], int(dimyx[1])),
               y=np.linspace(lat_range[0], lat_range[1], int(dimyx[0])),
               _x=gx, _y=gy,
               bkgd=out['bkgd'], total=out['total'],
               clust=out['clust'], lamb=out['lamb'])
    return out



# output rates algorithm (original c++)
def cxxrates(param, revents, bwd, tperiod, gx, gy):
    t = np.array(revents['tt'])
    x = np.array(revents['xx'])
    y = np.array(revents['yy'])
    m = np.array(revents['mm'])
    pb = np.array(revents['prob'])
    bwd = np.array(bwd)
    
    # extract time period information
    tstart2 = tperiod['study_start']
    tlength = tperiod['study_end']
    
    mu = param['mu']
    A = param['A']
    c = param['c']
    alpha = param['alpha']
    p = param['p']
    D = param['D']
    q = param['q']
    gamma = param['gamma']
    
    N = len(t)
    ngx = gx.shape[0]
    ngy = gy.shape[0]
    
    bkgd = np.zeros([ngx, ngy])
    total = np.zeros([ngx, ngy])
    clust = np.zeros([ngx, ngy])
    lamb = np.zeros([ngx, ngy])

    for i in range(0, ngx):
        for j in range(0, ngy):
            tmp = dGauss(dist(x, y, gx[i], gy[j]), bwd)
            sum1 = np.sum(pb * tmp)
            sum2 = np.sum(tmp)
            # # old loop
            # sum1 = 0.
            # sum2 = 0.
            # for l in range(0, N):
            #     tmp = dGauss(dist(x[l], y[l], gx[i], gy[j]), bwd[l])
            #     sum1 = sum1 + pb[l] * tmp
            #     sum2 = sum2 + tmp
            bkgd[i, j] = sum1 / (tlength - tstart2)
            total[i, j] = sum2 / (tlength - tstart2)
            clust[i, j] = 1. - sum1 / sum2
            lamb[i, j] = mu * bkgd[i, j]
            
            tmp2 = A * np.exp(alpha * m) * (p - 1)/c * \
                    np.power(1 + (tlength - t)/c, - p) * \
                    (q - 1) / (D * np.exp(gamma * m) * np.pi) * \
                    np.power(1 + dist2(x, y, gx[i], gy[j]) / (D * np.exp(gamma * m)), - q)
            lamb[i, j] = lamb[i, j] + np.sum(tmp2)
            # # old loop
            # for l in range(0, N):
            #     lamb[i, j] += A * np.exp(alpha * m[l]) * \
            #                   (p - 1)/c * pow(1 + (tlength - t[l])/c, - p) * \
            #                   (q - 1) / (D * np.exp(gamma * m[l]) * np.pi) * \
            #                   pow(1 + dist2(x[l], y[l], gx[i], gy[j]) / (D * np.exp(gamma * m[l])), - q)
    return dict(bkgd = bkgd, total = total, clust = clust, lamb = lamb)



def plot_rates(rts, fit):
    
    figsize = (6,6)

    extent = [min(rts['x']), max(rts['x']),
              min(rts['y']), max(rts['y']),]

    titles = ['Background seismicity rate',
              'Total spatial seismicity rate',
              'Clustering coefficient',
              'Conditional intensity function (end of study period)']
    
    keys = list(rts.keys())[-4:]
    
    yy, xx = np.meshgrid(rts['y'], rts['x'])
    
    lon_step = 2
    lat_step = 2
    
    fig = plt.figure(figsize=figsize)
    
    for i in range(0,4):
        
        ax = plt.subplot(2,2,i+1)
        
        bm = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[2],
                     urcrnrlon=extent[1], urcrnrlat=extent[3],
                     projection='cyl', resolution='l', fix_aspect=False, ax=ax)
        bm.drawcoastlines()
        bm.drawcountries()
        
        data = rts[keys[i]]
        colormesh = bm.pcolormesh(xx, yy, data, cmap='jet') # contour
        bm.colorbar(colormesh) #, label=titles[i])     
        ax.set_title(titles[i])
        
        shape = [fit.catalog.region_poly['lon'], fit.catalog.region_poly['lat']]
        patches = [Polygon(np.array(shape).transpose(), True, edgecolor='r',
                           fill=False, linewidth=1.5, zorder=20)]
        ax.add_collection(PatchCollection(patches, match_original=True))

        if i >= 2:
            ax.set_xlabel('Longitude (°)') #, labelpad=15)
        if i % 2 == 0:
            ax.set_ylabel('Latitude (°)') #, labelpad=25)                
        
        # Add meridian and parallel gridlines
        meridians = np.round(np.arange(extent[0], extent[1] + lon_step, lon_step), 2)
        parallels = np.round(np.arange(extent[2], extent[3] + lat_step, lat_step), 2)
        ax.set_yticks(parallels)
        ax.set_yticklabels(parallels, verticalalignment='center', horizontalalignment='right')
        ax.set_xticks(meridians)
        ax.set_xticklabels(meridians)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.grid(False)
            
    plt.show()



def plot_map_events(prob, fit):
    # plot triggered and background events 95% confidence
    
    df_all = fit.catalog.revents
    proj = CatalogueEtas.xy2longlat(df_all['xx'], df_all['yy'],
                                    region_poly=fit.catalog.region_poly,
                                    dist_unit=fit.catalog.dist_unit)
    df_all['longitude'] = proj['long']
    df_all['latitude'] = proj['lat']
    
    df_triggered = df_all[ (df_all['flag']==1) & (prob['prob']>0.95) ]
    df_background = df_all[ (df_all['flag']==1) & ((1-prob['prob'])>0.95) ]
    
    region = fit.catalog.region_win
    region_poly = fit.catalog.region_poly
    
    extent = [164, 181, -48, -32]
    lon_step = 1
    lat_step = 1
    
    figsize = (6,6)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    ax.scatter(df_all['longitude'], df_all['latitude'], (4+df_all['mm']),
                color=[0.7,0.7,0.7], alpha=0.5, label='Events until 2010')
    ax.scatter(df_triggered['longitude'], df_triggered['latitude'], (4+df_triggered['mm']),
                color='r', alpha=0.5, label='Triggered events with 95% confidence')
    ax.scatter(df_background['longitude'], df_background['latitude'], (4+df_background['mm']),
                color='b', alpha=0.5, label='Background events with 95% confidence')

    bm = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[2],
                  urcrnrlon=extent[1], urcrnrlat=extent[3],
                  projection='cyl', resolution='l', fix_aspect=False, ax=ax)
    bm.drawcoastlines()
    bm.drawcountries()
    
    shape = [region_poly['lon'], region_poly['lat']]
    from matplotlib.patches import Polygon
    patches = [Polygon(np.array(shape).transpose(), True, edgecolor='r',
                        fill=False, linewidth=1.5, zorder=20)]
    ax.add_collection(PatchCollection(patches, match_original=True))    
    
    
    ax.set_xlabel('Longitude (°)') #, labelpad=15)
    ax.set_ylabel('Latitude (°)') #, labelpad=25)
    
    # Add meridian and parallel gridlines
    meridians = np.round(np.arange(extent[0], extent[1] + lon_step, lon_step), 2)
    parallels = np.round(np.arange(extent[2], extent[3] + lat_step, lat_step), 2)
    ax.set_yticks(parallels)
    ax.set_yticklabels(parallels, verticalalignment='center', horizontalalignment='right')
    ax.set_xticks(meridians)
    ax.set_xticklabels(meridians)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.grid(False)
    plt.legend()
        
    plt.show()
    
        
    
    


#%%

if __name__ == "__main__":
    
    etas_iran = load_pickle('test/etas_iran')
    
    pr = probs(etas_iran)
    pr['prob'].describe()
    
    r = rates(etas_iran)
    