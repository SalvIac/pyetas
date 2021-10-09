# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import numpy as np
import pandas as pd
from openquake.hmtk.seismicity.declusterer.dec_gardner_knopoff import GardnerKnopoffType1
from openquake.hmtk.seismicity.declusterer.distance_time_windows import GardnerKnopoffWindow
from pyetas.gardner_knopoff_window import GardnerKnopoffWindowOrig
from pyetas.gardner_knopoff import GardnerKnopoff


def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a


def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])


def find_ellipse(x, y):
    xmean = x.mean()
    ymean = y.mean()
    x -= xmean
    y -= ymean
    a = fitEllipse(x,y)
    center = ellipse_center(a)
    center[0] += xmean
    center[1] += ymean
    phi = ellipse_angle_of_rotation(a)
    axes = ellipse_axis_length(a)
    x += xmean
    y += ymean
    return center, phi, axes


def decluster_ss(catalogue):
    config = {'time_distance_window': GardnerKnopoffWindow(),
              'fs_time_prop': 1.0}   
    gk = GardnerKnopoffType1()
    vcl, flagvector = gk.decluster(catalogue, config)
    return vcl, flagvector


def decluster(catalogue):
    config = {'time_distance_window': GardnerKnopoffWindowOrig(),
              'fs_time_prop': 1.0}   
    gk = GardnerKnopoff()
    vcl, flagvector = gk.decluster(catalogue, config)
    return vcl, flagvector


def get_cluster_info(catalogue, vcl, flagvector):
    # count cluster events
    num_events = np.array([np.sum(vcl==i) for i in range(0, np.max(vcl))])
    cluster = np.array(range(0, np.max(vcl)))
    # magnitude of mainshock cluster
    mags = list()
    datetimes = list()
    lats = list()
    lons = list()
    num_aftershocks = list()
    publicid = list()
    for i in range(0, np.max(vcl)):
        ind = np.logical_and(vcl==i, flagvector==0)
        if np.sum(ind) == 1:
            mags.append(catalogue.data['magnitude'][ind][0])
            datetimes.append((catalogue.data['datetime'][ind][0]))
            lons.append(catalogue.data['longitude'][ind][0])
            lats.append(catalogue.data['latitude'][ind][0])
            num_aftershocks.append(np.sum((vcl==cluster[i]) &
                                          (catalogue.data['datetime'] > datetimes[-1])))
            publicid.append(catalogue.data['publicid'][ind][0])
        else:
            mags.append(np.nan)
            datetimes.append(np.datetime64("1000-01-01"))
            lons.append(np.nan)
            lats.append(np.nan)
            num_aftershocks.append(np.nan)
            publicid.append(np.nan)
    return pd.DataFrame(zip(cluster, num_events, num_aftershocks, mags, lons,
                            lats, datetimes, publicid),
                        columns=['cluster', 'num_events', "num_aftershocks",
                                 'magnitude', "longitude", "latitude",
                                 "datetime", "publicid"])


def get_default_bkg(lonsx, latsy, bin=0.1):
    lonsx = np.arange(np.floor(lonsx[0]*10)/10-bin,
                      np.ceil(lonsx[1]*10)/10+bin+0.1, bin)
    latsy = np.arange(np.floor(latsy[0]*10)/10-bin,
                      np.ceil(latsy[1]*10)/10+bin+0.1, bin)
    rated = np.ones((lonsx.shape[0], latsy.shape[0]))
    bkg = {"type": "gridded",
           'lon': lonsx, 'lat': latsy,
           'rated': rated*0.1}
    return bkg

