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

import numpy as np

from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point
from shapely.errors import TopologicalError


    
def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.
    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.
    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)




def get_voronoi(rpoly, min_prec=0.05):
    try:
        coords = [(x, y) for x,y in zip(rpoly["px"], rpoly["py"])]
        region_win = Polygon(coords)
        
        # mesh of points
        precision = np.min([min_prec, region_win.length/1000])
        print("precision Voronoi:", precision)
        x = np.arange( np.floor(np.min(rpoly["px"])/precision)*precision,
                       np.ceil(np.max(rpoly["px"])/precision)*precision,
                       precision)
        y = np.arange( np.floor(np.min(rpoly["py"])/precision)*precision,
                       np.ceil(np.max(rpoly["py"])/precision)*precision,
                       precision)    
        xv, yv = np.meshgrid(x, y)
        points = np.vstack([xv.flatten(), yv.flatten()]).transpose()
        
        # points inside larger polygon
        flag = np.array([Point(xxx, yyy).within(region_win) for xxx, yyy in
                          zip(points[:,0], points[:,1])])
        points = points[flag,:]
        # points_orig = points.copy()
        
        # iterate twice to change the points location to the centroids (more precise)
        for _ in range(0,2):
            vor = Voronoi(points) # compute Voronoi tesselation
            regions, vertices = voronoi_finite_polygons_2d(vor) # get regions and vertices
            regions_vert = []
            areas = []
            points = []
            for region in regions:
                polygon = vertices[region]
                # clipping polygon
                poly = Polygon(polygon)
                poly = poly.intersection(region_win)
                # area, centroids and regions
                areas.append( poly.area )
                points.append( [poly.centroid.xy[0][0], poly.centroid.xy[1][0]] )
                regions_vert.append( [p for p in poly.exterior.coords] )
            points = np.array(points)
        areas = np.array(areas)
    except TopologicalError:
        print("min_prec halved!")
        return get_voronoi(rpoly, min_prec/2)
    return points, areas, regions_vert





def get_mixed_voronoi(rpoly, x0, min_prec=0.05, disc_num=100, shape_factor=2):
    try:
        coords = [(x, y) for x,y in zip(rpoly["px"], rpoly["py"])]
        region_win = Polygon(coords)
        
        # mesh of points
        xx, yy = region_win.exterior.coords.xy
        dist = np.array([Point(xxx, yyy).distance(Point(x0[0], x0[1]))
                         for xxx, yyy in zip(xx, yy)])
        max_dist = np.ceil(np.max(dist))
        
        r = np.logspace(-10, np.log10(max_dist), disc_num)
        # r2 = np.linspace(0., max_dist, disc_num)
        
        # ind = np.diff(r1) < np.diff(r2)
        # ind = np.concatenate( ([True], ind) )
        # dif = np.concatenate( (np.diff(r1)[ind], np.diff(r2)[~ind]) )
        # r = np.concatenate( (r1[ind], r2[~ind]) )
        theta = np.linspace(0, 2*np.pi, int(disc_num/shape_factor))
        rv, thetav = np.meshgrid(r, theta)
        
        
        x = x0[0] + rv.flatten()*np.cos(thetav.flatten())
        y = x0[1] + rv.flatten()*np.sin(thetav.flatten())
        points1 = np.vstack([x, y]).transpose()
        points1 = np.unique(points1, axis=0)
        app = np.array( [[x0[0] + max_dist+1, x0[1] + max_dist+1],
                         [x0[0] + max_dist+1, x0[1] - max_dist+1],
                         [x0[0] - max_dist+1, x0[1] + max_dist+1],
                         [x0[0] - max_dist+1, x0[1] - max_dist+1]] )
        points1 = np.concatenate( (points1, app), axis=0 )
        points1 = np.concatenate( (points1, np.array( [x0])), axis=0 )
        
        
        # mesh of points
        precision = min_prec # np.min([min_prec, region_win.length/disc_num])
        print("precision Voronoi:", precision)
        x = np.arange( np.floor(np.min(rpoly["px"])/precision)*precision,
                       np.ceil(np.max(rpoly["px"])/precision)*precision,
                       precision)
        y = np.arange( np.floor(np.min(rpoly["py"])/precision)*precision,
                       np.ceil(np.max(rpoly["py"])/precision)*precision,
                       precision)    
        xv, yv = np.meshgrid(x, y)
        points2 = np.vstack([xv.flatten(), yv.flatten()]).transpose()
        
        # points inside larger polygon
        flag = np.array([Point(xxx, yyy).within(region_win) for xxx, yyy in
                          zip(points2[:,0], points2[:,1])])
        points2 = points2[flag,:]

        points = np.vstack((points1, points2))


        # iterate twice to change the points location to the centroids (more precise)
        for i in range(0,2):
            vor = Voronoi(points) # compute Voronoi tesselation
            regions, vertices = voronoi_finite_polygons_2d(vor) # get regions and vertices
            regions_vert = []
            areas = []
            points = []
            for region in regions:
                polygon = vertices[region]
                # clipping polygon
                poly = Polygon(polygon)
                poly = poly.intersection(region_win)
                if poly.area != 0.:
                    # area, centroids and regions
                    areas.append( poly.area )
                    points.append( [poly.centroid.xy[0][0], poly.centroid.xy[1][0]] )
                    regions_vert.append( [p for p in poly.exterior.coords] )
            points = np.array(points)
            if i == 0:
                points = np.concatenate( (points, app), axis=0 )
        areas = np.array(areas)
    except TopologicalError:
        print("disc_num doubled!")
        return get_mixed_voronoi(rpoly, x0, min_prec/2, disc_num*2)
    return points, areas, regions_vert






def get_polar_voronoi(rpoly, x0, disc_num=100, shape_factor=2):
    try:
        coords = [(x, y) for x,y in zip(rpoly["px"], rpoly["py"])]
        region_win = Polygon(coords)
        
        # mesh of points
        xx, yy = region_win.exterior.coords.xy
        dist = np.array([Point(xxx, yyy).distance(Point(x0[0], x0[1]))
                         for xxx, yyy in zip(xx, yy)])
        max_dist = np.ceil(np.max(dist))
        
        r = np.logspace(-10, np.log10(max_dist), disc_num)
        # r2 = np.linspace(0., max_dist, disc_num)
        
        # ind = np.diff(r1) < np.diff(r2)
        # ind = np.concatenate( ([True], ind) )
        # dif = np.concatenate( (np.diff(r1)[ind], np.diff(r2)[~ind]) )
        # r = np.concatenate( (r1[ind], r2[~ind]) )
        theta = np.linspace(0, 2*np.pi, int(disc_num/shape_factor))
        rv, thetav = np.meshgrid(r, theta)
        
        
        x = x0[0] + rv.flatten()*np.cos(thetav.flatten())
        y = x0[1] + rv.flatten()*np.sin(thetav.flatten())
        points = np.vstack([x, y]).transpose()
        points = np.unique(points, axis=0)
        app = np.array( [[x0[0] + max_dist+1, x0[1] + max_dist+1],
                         [x0[0] + max_dist+1, x0[1] - max_dist+1],
                         [x0[0] - max_dist+1, x0[1] + max_dist+1],
                         [x0[0] - max_dist+1, x0[1] - max_dist+1]] )
        points = np.concatenate( (points, app), axis=0 )
        points = np.concatenate( (points, np.array( [x0])), axis=0 )
        


        
        # points inside larger polygon
        # flag = np.array([Point(xxx, yyy).within(region_win) for xxx, yyy in
        #                   zip(points[:,0], points[:,1])])
        # points = points[flag,:]
        # points_orig = points.copy()
        
        # iterate twice to change the points location to the centroids (more precise)
        for i in range(0,2):
            vor = Voronoi(points) # compute Voronoi tesselation
            regions, vertices = voronoi_finite_polygons_2d(vor) # get regions and vertices
            regions_vert = []
            areas = []
            points = []
            for region in regions:
                polygon = vertices[region]
                # clipping polygon
                poly = Polygon(polygon)
                poly = poly.intersection(region_win)
                if poly.area != 0.:
                    # area, centroids and regions
                    areas.append( poly.area )
                    points.append( [poly.centroid.xy[0][0], poly.centroid.xy[1][0]] )
                    regions_vert.append( [p for p in poly.exterior.coords] )
            points = np.array(points)
            if i == 0:
                points = np.concatenate( (points, app), axis=0 )
        areas = np.array(areas)
    except TopologicalError:
        print("disc_num doubled!")
        return get_polar_voronoi(rpoly, x0, disc_num*2)
    return points, areas, regions_vert






#%%


# if __name__ == "__main__":
    
#     import matplotlib.pyplot as plt
#     from myutils.utils_pickle import load_pickle, save_pickle
#     from scipy.spatial import voronoi_plot_2d
    
    
#     theta = load_pickle('test/param1')
#     rdata = load_pickle('test/rdata')
#     revents = rdata['revents']
#     rpoly = rdata['rpoly']
#     tperiod = rdata['tperiod']
#     integ0 = rdata['integ0']
#     ihess = load_pickle('test/ihess')
#     rverbose = verbose = 1    
    
    
#     # points, areas, regions_vert = get_voronoi(rpoly, min_prec=0.05)
#     # # plot
#     # fig = plt.figure()
#     # for polygon in regions_vert:
#     #     plt.fill(*zip(*polygon), alpha=0.4, color=np.random.random(3))
#     # plt.plot(points[:, 0], points[:, 1], 'ro')
#     # plt.axis('equal')
#     # plt.show()
    
#     # x0 = [-1., 2]
#     # points, areas, regions_vert = get_polar_voronoi(rpoly, x0, disc_num=100)
#     # # plot
#     # fig = plt.figure()
#     # for polygon in regions_vert:
#     #     plt.fill(*zip(*polygon), alpha=0.4, color=np.random.random(3))
#     # plt.plot(points[:, 0], points[:, 1], 'ro')
#     # plt.axis('equal')
#     # plt.show()
    
    
        
#     x0 = [-1., 2]
#     points, areas, regions_vert = get_mixed_voronoi(rpoly, x0, min_prec=1,
#                                                     disc_num=20)
#     # plot
#     fig = plt.figure()
#     for polygon in regions_vert:
#         plt.fill(*zip(*polygon), alpha=0.4, color=np.random.random(3))
#     plt.plot(points[:, 0], points[:, 1], 'ro')
#     plt.axis('equal')
#     plt.show()
    
    
    
