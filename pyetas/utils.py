# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import numpy as np
from shapely.geometry import Point, Polygon
from openquake.hmtk.seismicity.catalogue import Catalogue
from pyetas.etas8p.catalog import CatalogueEtas




def get_region(xvals, yvals):
    coords = [(x, y) for x,y in zip(xvals, yvals)]
    region_win = Polygon(coords)
    return region_win



def longlat2xy(long, lat, centroid_lon, centroid_lat, dist_unit="degree"):
    # Equirectangular projection
    if len(long) != len(lat):
        raise Exception("long and lat must be of equal length.")
    if dist_unit == "degree":
        x = np.cos(centroid_lat * np.pi / 180) * (long - centroid_lon)
        y = lat - centroid_lat
    elif dist_unit == "km":
        x = 111.320 * np.cos(lat/180 * np.pi) * long
        y = 110.574 * lat
    else:
        raise Exception("dist.unit argument must be either degree or km.")
    return dict(x=x, y=y)



def xy2longlat(x, y, centroid_lon, centroid_lat, dist_unit="degree"):
    if len(x) != len(y):
        raise Exception("x and y must be of equal length.")
    if dist_unit == "degree":
        long = x / np.cos(centroid_lat * np.pi / 180) + centroid_lon
        lat = y + centroid_lat
    elif dist_unit == "km":
        lat = y / 110.574
        long = x / (111.320 * np.cos(lat/180 * np.pi))
    else:
        raise Exception("dist.unit argument must be either degree or km.")
    return dict(long=long, lat=lat)




def project(lon, lat, fit):
    proj = CatalogueEtas.longlat2xy(long=np.array(lon),
                                    lat=np.array(lat),
                                    region_poly=fit.catalog.region_poly,
                                    dist_unit=fit.catalog.dist_unit)
    return proj['x'], proj['y']


def get_projected_region(lon, lat, fit):
    projx, projy = project(lon, lat, fit)
    return get_region(projx, projy)


def get_oq_catalogue(df):
    # openquake object
    data = dict()
    for key in df.columns:
        data[key] = df[key].to_numpy()
    return Catalogue.make_from_dict(data)


def filter_polygon(df, poly):
    region_win = get_region(poly[:,0], poly[:,1])
    flag = np.array([Point(xxx, yyy).within(region_win) for xxx, yyy in 
                     zip(df['longitude'], df['latitude'])])
    df = df[flag]
    return df


def get_buffer_region(region, dist=0.5):
    buf_region = {'lon': np.array([np.min(region['lon'])-dist, np.max(region['lon'])+dist,
                                    np.max(region['lon'])+dist, np.min(region['lon'])-dist]),
                  'lat': np.array([np.max(region['lat'])+dist, np.max(region['lat'])+dist,
                                    np.min(region['lat'])-dist, np.min(region['lat'])-dist])}
    return buf_region




