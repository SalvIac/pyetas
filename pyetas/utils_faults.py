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
"""

import re
import numpy as np
from openquake.hazardlib.geo.mesh import RectangularMesh, Mesh



def read_surface_fsp(filename):
    with open(filename, 'r') as reader:
        text = reader.readlines()
    lons = []
    lats = []
    depths = []
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



def read_metadata_fsp(filename):
    with open(filename, 'r') as reader:
        text = reader.readlines()
    line = text[2].split('\t')
    year = float(line[2].split("/")[2])
    month = float(line[2].split("/")[0])
    day = float(line[2].split("/")[1])
    mag = float( text[6].split("Mw = ")[1].split('\t')[0] )
    return dict(year=year, month=month, day=day, mag=mag)



def read_surface_fsp_mesh(filename):
    with open(filename, 'r') as reader:
        text = reader.readlines()
    lons = []
    lats = []
    depths = []
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
    index = []
    for line in text:
        try:
            new_list = list(filter(None, re.split(' |\t', line[:-1]))) # line[:-1].split(' ')
            fl_list = [float(el) for el in new_list]
            index.append(fl_list[-1])
        except:
            continue
    lons = []
    lats = []
    depths = []
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
    
    lons_2d = []
    lats_2d = []
    depths_2d = []
    for depth in np.unique(depths):
        lons_2d.append(lons[depths == depth])
        lats_2d.append(lats[depths == depth])
        depths_2d.append(depths[depths == depth])
    lons_2d = np.array(lons_2d, dtype=object)
    lats_2d = np.array(lats_2d, dtype=object)
    depths_2d = np.array(depths_2d, dtype=object)
    
    check = lons.shape[0] != lons_2d.shape[0]*lons_2d.shape[1]
    if check:
        raise Exception('error with lat lon depth')
    
    return RectangularMesh(lons_2d, lats_2d, depths_2d)

