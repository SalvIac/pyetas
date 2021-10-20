# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""


import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection



class MapEarthquakes():
    
    figsize = (6,6)
    lon_step = 2
    lat_step = 2
    resolution = 'l'


    def __init__(self, extent, ax=None, figsize=None):
        # Define lower left, uperright lontitude and latitude respectively
        self.extent = extent
        if figsize is not None:
            self.figsize = figsize
        if ax is None:
            self.fig = plt.figure(figsize=self.figsize)
            self.ax = self.fig.gca()
        else:
            self.ax = ax

        
    def plot_base(self):
        # Create a basemap instance that draws the Earth layer
        bm = Basemap(llcrnrlon=self.extent[0], llcrnrlat=self.extent[2],
                     urcrnrlon=self.extent[1], urcrnrlat=self.extent[3],
                     projection='cyl', resolution=self.resolution,
                     fix_aspect=False, ax=self.ax)
        bm.drawcoastlines()
        bm.drawcountries()
        self.bm = bm
    

    def add_events(self, df, **kwargs):
        self.ax.scatter(df['longitude'], df['latitude'], df['magnitude'] ** 2, 
                        **kwargs)


    def add_polygon(self, poly, **kwargs):
        # poly is a (n,2) numpy array
        patches = [Polygon(poly, True, **kwargs)]
        self.add_patches(patches)
        

    def add_patches(self, patches):
        if not isinstance(patches, list):
            patches = list(patches)
        self.ax.add_collection(PatchCollection(patches, match_original=True))


    def look(self):        
        self.ax.set_xlabel('Longitude (°)') #, labelpad=15)
        self.ax.set_ylabel('Latitude (°)') #, labelpad=25)
        # Add meridian and parallel gridlines
        meridians = np.round(np.arange(self.extent[0], self.extent[1]
                                       + self.lon_step, self.lon_step), 2)
        parallels = np.round(np.arange(self.extent[2], self.extent[3]
                                       + self.lat_step, self.lat_step), 2)
        self.ax.set_yticks(parallels)
        self.ax.set_yticklabels(parallels, verticalalignment='center',
                                horizontalalignment='right')
        self.ax.set_xticks(meridians)
        self.ax.set_xticklabels(meridians)
        self.ax.set_xlim(self.extent[0], self.extent[1])
        self.ax.set_ylim(self.extent[2], self.extent[3])
        self.ax.grid(False)
    
    
    def legend(self):
        self.ax.legend()
    
        
    def get_fig_ax(self):    
        return self.fig, self.ax
        
   
    def simple_save(self, filename=None):
        self.fig.tight_layout()    
        plt.savefig(filename, bbox_inches='tight')
        self.fig.clear()
        plt.close(self.fig)
        
    
    def save(self, filename=None, crop=False):
        filename = self.adjust_filename(filename)
        self.simple_save(filename)
        if crop:
            self.crop_image(filename)
    
    
    @classmethod
    def adjust_filename(cls, filename):
        if filename is None:
            return 'test.png'
        if '.' not in filename:
            return filename+'.png'
        return filename
            
            
    @classmethod
    def crop_image(cls, filename, coeff=None):
        if coeff is None:
            coeff = cls.coeff
        # crop image
        im = Image.open(filename)
        width, height = im.size
        x = coeff[0]*width
        y = coeff[1]*height
        widthx = coeff[2]*width
        heighty = coeff[3]*height
        cropped_image = im.crop((x, y, widthx, heighty))
        cropped_image.save(filename)
    
    
    def show(self):
        self.fig.tight_layout()    
        plt.show()
    
