# -*- coding: utf-8 -*-
'''
Created on Thu May  5 08:51:36 2022

@author: alhe551

Serie of mapping function used during PhD to ensure homogemeous display.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
os.environ['PROJ_LIB'] = os.path.join(os.environ['CONDA_PREFIX'], 'share', 'proj')
from mpl_toolkits.basemap import Basemap


#%% class zone
class zone(object):
    def __init__(self,area):
        '''
        Selection of the area for mapping
        Input:
            area:            str() from list: HG, HG-BoP, BckBn
        Output:
            LOWER_LEFT_LON:  int()
            LOWER_LEFT_LAT   int()
            UPPER_RIGHT_LON: int()
            UPPER_RIGHT_LAT: int()
        WARNING: unknown/default location 10*10deg square around point nemo.
        '''
        self.area = area
        if area == 'HG-BoP':
            self.LOWER_LEFT_LAT = -38.033157419989735
            self.LOWER_LEFT_LON = 174.54998618909264
            self.UPPER_RIGHT_LAT = -35.993110976936805
            self.UPPER_RIGHT_LON = 178.7723699238572
        elif area == 'HG':
            self.LOWER_LEFT_LAT = -37.24822666266156
            self.LOWER_LEFT_LON = 174.4191057667842
            self.UPPER_RIGHT_LAT = -35.664520480112174
            self.UPPER_RIGHT_LON = 176.18599269342826
        elif area == 'BckBn':
            self.LOWER_LEFT_LAT = -51.981240790604566
            self.LOWER_LEFT_LON = 161.03022670026627
            self.UPPER_RIGHT_LAT = -31.026108847126423
            self.UPPER_RIGHT_LON = 184.96977329971511
        else:
            self.area = 'Point Nemo'
            self.LOWER_LEFT_LAT = -43.876667
            self.LOWER_LEFT_LON = -128.393333
            self.UPPER_RIGHT_LAT = -52.876667
            self.UPPER_RIGHT_LON = -118.393333

#%% corner_remind
    def corner_remind(self):
        return print(f'{self.area} area coordinate:\n{self.LOWER_LEFT_LON}, {self.LOWER_LEFT_LAT},\n{self.UPPER_RIGHT_LON}, {self.UPPER_RIGHT_LAT}')

#%% cut_corners
    def cut_corners(self, data_nc):
        '''
        Cut the netcdf dataset to the dimension set by map_corners
        Input:
            data_nc: .nc dataset/dataArray from CMEMS
        Output:
            new_nc:  dataset/dataArray restricted by long/lat set    
        '''
        new_nc = data_nc.sel(lat=slice(self.UPPER_RIGHT_LAT,self.LOWER_LEFT_LAT)).sel(lon=slice(self.LOWER_LEFT_LON, self.UPPER_RIGHT_LON))
        return new_nc

#%% map_raster
    def map_raster(self, data_nc, min_col=-1e5, max_col=-1e5,
               legend = 'legend',
               map_col = 'viridis'):
        '''
        Create and save a map using the data provided for log de CHL
        Input:
            data_nc:    .nc file from copernicus (CMEMS)
            min_col:    minimal value for basemap.colormesh
            max_col:    maximal value for basemap.colormesh
            legend:     str() object
            map_col:    cmap name in str()
        Output:
            fig:        figure object
        '''
        if min_col == -1e5 and max_col == -1e5:
            min_col = np.nanmin(data_nc)
            max_col = np.nanmax(data_nc)
            
        plt.figure(figsize=(12,8)) #figsize=(length, height)
        fig, axs = plt.subplots() #(row, col)

        mp = Basemap(projection='merc', # projection.
                     #List available: https://yonsci.github.io/yon_academic/portfolio/portfolio-8/
                     llcrnrlon = self.LOWER_LEFT_LON,    # lower longitude
                     llcrnrlat = self.LOWER_LEFT_LAT,   # lower latitude
                     urcrnrlon = self.UPPER_RIGHT_LON,   # uppper longitude
                     urcrnrlat = self.UPPER_RIGHT_LAT,  # uppper latitude
                     resolution = 'f')
        # resolutions: c - crude, l - low, i - intermediate, h - high, f - full

        x,y = mp(*np.meshgrid(data_nc.lon, data_nc.lat)) # longitude-180 if over 360°
        var_raster = mp.pcolormesh(x,y, data_nc, vmin=min_col, vmax=max_col,
                                   cmap=map_col)
        cb = mp.colorbar(var_raster,'right', size='3%', pad='2%')#size=thickness, pad=distance to figure
        cb.set_label(legend, fontsize=9)

        mp.drawcoastlines(linewidth=.5) # draw coast based on mp
        mp.fillcontinents(color='grey') # draw continent based on mp
        mp.drawmapboundary(fill_color='white') # fill the background of the map
        #mp.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # pixel style drawing ?

        parallels = np.arange(int(self.LOWER_LEFT_LAT)-1,int(self.UPPER_RIGHT_LAT)+1,.5)
        meridians = np.arange(int(self.LOWER_LEFT_LON)-1,int(self.UPPER_RIGHT_LON)+1,.5)
        mp.drawparallels(parallels,labels=[1,0,0,0],fontsize=10) #labels [left, right, ?, ?]
        mp.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10) #labels [?, ?, up, down]

        #axs.set_xlabel('Test') #resolve the overlap issue
        #axs.set_ylabel('Test2') #resolve the overlap issue
        #axs.set_title(figure_title)
        fig.tight_layout()

        plt.close(fig)
        
        return fig

#%% map_contour
    def map_contour(self, data_nc, min_col=-1e5, max_col=-1e5,
               legend = 'legend',
               lin_col = 'black'):
        '''
        Create and save a map using the data provided for log de CHL
        Input:
            data_nc:    .nc file from copernicus (CMEMS)
            min_col:    minimal value for basemap.colormesh
            max_col:    maximal value for basemap.colormesh
            legend:     str() object
            map_col:    cmap name in str()
        Output:
            fig:        figure object
        '''
        if min_col == -1e5 and max_col == -1e5:
            min_col = np.nanmin(data_nc)
            max_col = np.nanmax(data_nc)

        plt.figure(figsize=(12,8)) #figsize=(length, height)
        fig, axs = plt.subplots() #(row, col)

        mp = Basemap(projection='merc', # projection.
                     #List available: https://yonsci.github.io/yon_academic/portfolio/portfolio-8/
                     llcrnrlon = self.LOWER_LEFT_LON,    # lower longitude
                     llcrnrlat = self.LOWER_LEFT_LAT,   # lower latitude
                     urcrnrlon = self.UPPER_RIGHT_LON,   # uppper longitude
                     urcrnrlat = self.UPPER_RIGHT_LAT,  # uppper latitude
                     resolution = 'f')
        # resolutions: c - crude, l - low, i - intermediate, h - high, f - full

        x,y = mp(*np.meshgrid(data_nc.lon, data_nc.lat)) # longitude-180 if over 360°
        var_raster = mp.contour(x, y, data_nc, vmin=min_col, vmax=max_col,
                                 #levels=np.linspace(min_col, max_col, 4),
                                 levels=(min_col, -.1, 0, .1, .2, max_col), colors=lin_col)
        axs.clabel(var_raster, var_raster.levels, inline=True, fontsize=10)

        mp.drawcoastlines(linewidth=.5) # draw coast based on mp
        mp.fillcontinents(color='grey') # draw continent based on mp
        mp.drawmapboundary(fill_color='white') # fill the background of the map
        #mp.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # pixel style drawing ?

        parallels = np.arange(int(self.LOWER_LEFT_LAT)-1,int(self.UPPER_RIGHT_LAT)+1,.5)
        meridians = np.arange(int(self.LOWER_LEFT_LON)-1,int(self.UPPER_RIGHT_LON)+1,.5)
        mp.drawparallels(parallels,labels=[1,0,0,0],fontsize=10) #labels [left, right, ?, ?]
        mp.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10) #labels [?, ?, up, down]

        #axs.set_xlabel('Test') #resolve the overlap issue
        #axs.set_ylabel('Test2') #resolve the overlap issue
        #axs.set_title(figure_title)
        fig.tight_layout()

        plt.close(fig)
        
        return fig

