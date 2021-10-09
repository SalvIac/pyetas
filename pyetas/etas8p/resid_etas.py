# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""


import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
plt.ioff()
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PatchCollection
from shapely.geometry import Point, Polygon

from pyetas.etas8p.etas import Etas
from pyetas.etas8p.from_mp import lambdaspatial, lambdatemporal, timetransform
from myutils.utils_pickle import load_pickle, save_pickle



def xy2longlat(x, y, region_poly, dist_unit="degree"):
    if len(x) != len(y):
        raise Exception("x and y must be of equal length.")
    region_poly_lon = np.array(region_poly['lon'])
    region_poly_lat = np.array(region_poly['lat'])
    if dist_unit == "degree":
        coords = [(x, y) for x,y in zip(region_poly_lon, region_poly_lat)]
        region_win = Polygon(coords)
        region_cnt = region_win.centroid.xy
        long = x / np.cos(region_cnt[1][0] * np.pi / 180) + region_cnt[0][0]
        lat = y + region_cnt[1][0]
    elif dist_unit == "km":
      lat = y / 110.574
      long = x / (111.320 * np.cos(lat/180 * np.pi) )
    else:
        raise Exception("dist.unit argument must be either degree or km.")
    return dict(long=long, lat=lat)




def resid_etas(fit, tpe="raw", n_temp=1000, dimyx=None):
    flg = np.array(fit.revents["flag"])
    tt = np.array(fit.revents["tt"])[flg == 1]
    xx = np.array(fit.revents["xx"])[flg == 1]
    yy = np.array(fit.revents["yy"])[flg == 1]
   
    tau = timetransform(fit)[flg == 1] - lambdatemporal([fit.catalog.rtperiod['study_start']], fit)
    U = 1 - np.exp(-np.diff(tau))
    tg = np.linspace(min(tt), max(tt), n_temp)
    tlam = lambdatemporal(tg, fit)

    # tt = [1,2,3,4,5,6,7,8,9,10]
    # tg = [4,2,7,100,10,25]

    # tt <- c(1,2,3,4,5,6,7,8,9,10)
    # tg <- c(2,4,7,9,25,100)
    # dfun = function(i){ sum((tt <= tg[i]) & (tt > tg[i - 1])) }
    
    # unlist(lapply(2:6, dfun)) 

    # tt = np.array([1,2,3,4,5,6,7,8,9,10])
    # tg = np.array([2,4,7,9,25,100])
    # np.array([sum((tt <= tg[i]) & (tt > tg[i - 1])) for i in range(1,6)])
    
    # unlist(lapply(2:n.temp, dfun))
    
    if tpe == 'raw':
        tres = np.array([sum((tt <= tg[i]) & (tt > tg[i - 1])) for i in range(1,n_temp)]) - tlam[1:] * np.diff(tg)
    elif tpe == 'reciprocal':
        tres = 1/np.tlam[1:] - np.diff(tg)
    elif tpe == 'pearson':
        tres = 1/np.sqrt(tlam[1:]) - np.sqrt(tlam[1:]) * np.diff(tg)
    
    # W = fit.catalog.region_win.exterior.xy
    # win = pointpats.Window( [(W[0][i], W[1][i]) for i in range(0, len(W[0]))] )
    
    W = fit.catalog.region_win
    xmin,ymin,xmax,ymax = W.bounds

    rv = np.diff([np.floor(xmin), np.ceil(xmax)])[0] / np.diff([np.floor(ymin), np.ceil(ymax)])[0]
    if (rv > 1):
        wide = int(round(1000.*rv))
        length = 1000
    else:
        wide = int(round(1000.*(1 / rv)))
        length = 1000

    xxx = np.linspace(np.floor(xmin), np.ceil(xmax), wide)
    yyy = np.linspace(np.floor(ymin), np.ceil(ymax), length)
    area = np.diff(xxx).min() * np.diff(yyy).min()
    
    # full grid
    xg_full, yg_full = np.meshgrid(xxx, yyy)
    xgg = xg_full.flatten()
    ygg = yg_full.flatten()
    # only points within polygon
    flag = np.array([(Point(xxx, yyy).within(W)) for xxx, yyy in zip(xgg,ygg)])
    xg = xgg[flag]
    yg = ygg[flag]
    wg = area*np.ones(yg.shape[0])
    proj = xy2longlat(xgg, ygg, region_poly=fit.catalog.region_poly,
                      dist_unit=fit.catalog.dist_unit)
    
    # fig = plt.figure()
    # ax = fig.add_subplot()
    # x,y = W.exterior.xy
    # ax.plot(x,y)
    # ax.scatter(xg, yg, s=1)
    # # ax.scatter(xg[flag], yg[flag], s=2)
    # plt.show()
    
    # Xs = spatstat::ppp(xx, yy, window=W, check=FALSE)
    # qd = spatstat::quadscheme(Xs)
    # xg = spatstat::x.quad(qd)
    # yg = spatstat::y.quad(qd)
    # wg = spatstat::w.quad(qd)
    # zg = spatstat::is.data(qd)
    
    # Xg = spatstat::ppp(xg, yg, window=W, check=FALSE)
    # spatstat::marks(Xg) = sres
    # sres = spatstat::Smooth(Xg, dimyx=dimyx, sigma=mean(fit$bwd))
    # gr = expand.grid(x=sres$xcol, y=sres$yrow)
    #  proj = xy2longlat(gr$x, gr$y, region.poly=fit$object$region.poly,
    #                    dist.unit=fit$object$dist.unit)
    # sres = data.frame(x=proj$long, y=proj$lat, z=c(t(sres$v)))
    # #sres = stats::na.omit(sres)

    # slam = lambdaspatial(xg, yg, fit)
    # zg = np.zeros(slam.shape[0])
    
    # if tpe == 'raw':
    #     _sres = zg - slam * wg
    # elif tpe == 'reciprocal':
    #     _sres = zg/slam - wg
    # elif tpe == 'pearson':
    #     _sres = zg/np.sqrt(slam) - np.sqrt(slam) * wg
        
    # distionary for sres
    sres=dict()
    # sres['x'] = np.empty(xg_full.shape)
    # sres['x'][:] = np.NaN
    # sres['y'] = np.empty(xg_full.shape)
    # sres['y'][:] = np.NaN
    # sres['z'] = np.empty(xg_full.shape)
    # sres['z'][:] = np.NaN
    # for i in range(0, xgg.shape[0]):
    #     to = np.where(np.logical_and(xg_full == xgg[i], yg_full == ygg[i]))
    #     sres['x'][to[0][0], to[1][0]] = proj['long'][i]
    #     sres['y'][to[0][0], to[1][0]] = proj['lat'][i]
    #     if flag[i]:
    #         fr = np.where(np.logical_and(xg == xgg[i], yg == ygg[i]))[0][0]
    #         to = np.where(np.logical_and(xg_full == xgg[i], yg_full == ygg[i]))
    #         sres['z'][to[0][0], to[1][0]] = _sres[fr]
    
    # if dimyx is None:
    #     rangex = [np.min(xg), np.max(xg)]
    #     rangey = [np.min(yg), np.max(yg)]
    #     rv = np.diff(rangex)[0] / np.diff(rangey)[0]
    #     if (rv > 1):
    #         dimyx = [128., round(128.*rv)]
    #     else:
    #         dimyx = [round(128.*(1 / rv)), 128.]

   
    out = dict(tau=tau, U=U, tres=tres, sres=sres, tpe=tpe, tg=tg)
    return out



def plot_resid(resid):

    U = resid['U']
    tau = resid['tau']
    tpe = resid['tpe']
    tres = resid['tres']
    tg = resid['tg']

    figsize = (6,6)
    fig = plt.figure(figsize=figsize)

    ax = plt.subplot(2,2,1)
    ax.plot(tg[1:], tres)
    ax.set_title(tpe+" temporal residuals")
    ax.plot([tg[1], tg[-1]], [0, 0], color='r', linestyle='--', linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel("Residuals")
    
    ax = plt.subplot(2,2,2)
    # ax.set_title(tpe+" spatial residuals")
    # extent = [np.min(sres['x']), np.max(sres['x']),
    #           np.min(sres['y']), np.max(sres['y'])]
    
    # bm = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[2],
    #               urcrnrlon=extent[1], urcrnrlat=extent[3],
    #               projection='cyl', resolution='h', fix_aspect=False, ax=ax)
    # bm.drawcoastlines()
    # bm.drawcountries()
    
    # zmax = np.nanmax(np.abs(sres['z']))
    # colormesh = bm.pcolormesh(sres['x'], sres['y'], sres['z'], cmap='jet') # contour
    # colormesh.set_clim(vmin=-0.43, vmax=0.43)
    # bm.colorbar(colormesh, label="spatial residuals")    
    
    # zmax = max(abs(sres['z']))
    # colormesh = bm.pcolormesh(sres['x'], sres['y'], sres['z'], cmap='jet') # contour
    # bm.colorbar(colormesh, vmin=-zmax, vmax=zmax, label="spatial residuals")     
    
    ax = plt.subplot(2,2,3)
    ax.plot(range(0, tau.shape[0]), tau, color='k')
    ax.plot(range(0, tau.shape[0]), range(0, tau.shape[0]), color='r')
    ax.set_title("transformed times")
    ax.set_xlabel("i")
    ax.set_ylabel("tau_i")
    ax.grid(True)
    
    ax = plt.subplot(2,2,4)
    stats.probplot(U, dist="uniform", plot=ax)
    ax.set_title("transformed times")
    ax.set_xlabel("observed~quantiles")
    ax.set_ylabel("quantiles~of~italic(U)(0,1)")
    ax.grid(True)

    plt.show()
   




#%%

if __name__ == "__main__":

    etas_iran = load_pickle('test/etas_iran')
    fit = etas_iran

    import time    
    
    res = resid_etas(etas_iran)
    plot_resid(res)
    # out = load_pickle('out')
    
