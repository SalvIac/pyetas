# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import re
import numpy as np
from openquake.hazardlib.geo.mesh import RectangularMesh, Mesh



def read_surface_fsp(filename):
    with open(filename, 'r') as reader:
        text = reader.readlines()
    lons = list()
    lats = list()
    depths = list()
    for line in text:
        if line[0] != '%':
            try:
                new_list = list(filter(None, line[:-1].split(' ')))
                fl_list = [float(el) for el in new_list]        
                lons.append(fl_list[1])
                lats.append(fl_list[0])
                depths.append(fl_list[4])
            except:
                continue
    return coords2recmesh(lons, lats, depths)



def read_surface_fsp_mesh(filename):
    with open(filename, 'r') as reader:
        text = reader.readlines()
    lons = list()
    lats = list()
    depths = list()
    for line in text:
        if line[0] != '%':
            try:
                new_list = list(filter(None, line[:-1].split(' ')))
                fl_list = [float(el) for el in new_list]        
                lons.append(fl_list[1])
                lats.append(fl_list[0])
                depths.append(fl_list[4])
            except:
                continue
    return Mesh(np.array(lons), np.array(lats), np.array(depths))



def read_surface_dat(filename):
    with open(filename, 'r') as reader:
        text = reader.readlines()
    index = list()
    for line in text:
        try:
            new_list = list(filter(None, re.split(' |\t', line[:-1]))) # line[:-1].split(' ')
            fl_list = [float(el) for el in new_list]
            index.append(fl_list[-1])
        except:
            continue
    lons = list()
    lats = list()
    depths = list()
    for line in text:
        try:
            new_list = list(filter(None, re.split(' |\t', line[:-1]))) # line[:-1].split(' ')
            fl_list = [float(el) for el in new_list]
            lons.append(fl_list[0])
            lats.append(fl_list[1])
            depths.append(fl_list[7])
        except:
            continue
    return Mesh(np.array(lons), np.array(lats), np.array(depths))

    


def coords2recmesh(lons, lats, depths):
    lons = np.array(lons)
    lats = np.array(lats)
    depths = np.array(depths)
    
    lons_2d = list()
    lats_2d = list()
    depths_2d = list()
    for depth in np.unique(depths):
        lons_2d.append(lons[depths == depth])
        lats_2d.append(lats[depths == depth])
        depths_2d.append(depths[depths == depth])
    lons_2d = np.array(lons_2d)
    lats_2d = np.array(lats_2d)
    depths_2d = np.array(depths_2d)
    
    check = lons.shape[0] != lons_2d.shape[0]*lons_2d.shape[1]
    if check:
        raise Exception('error with lat lon depth')
    
    return RectangularMesh(lons_2d, lats_2d, depths_2d)

