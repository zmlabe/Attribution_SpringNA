"""
Plot trends in FMA from 1979 for SST

Author    : Zachary M. Labe
Date      : 2 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'SST'
sliceperiod = 'April'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,lev1,data = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,'surface')

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2020
trend = UT.linearTrendR(np.asarray(data),years,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################
### Plot map of trends
limit = np.arange(-1,1.01,0.02)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{Trend [$^{\circ}$/decade]}'

ax = plt.subplot(111)

var = trend

m = Basemap(projection='ortho',lon_0=210,lat_0=0,resolution='l',area_thresh=1000)
m.drawcoastlines(color='darkgrey',linewidth=1)
m.drawstates(color='darkgrey',linewidth=0.5)
m.drawcountries(color='darkgrey',linewidth=0.5)
    
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,var,limit,extend='both')

cs1.set_cmap(cmocean.cm.balance)

m.fillcontinents(color='dimgrey')

plt.title(r'\textbf{%s %s-%s}' % (sliceperiod,yearmin,yearmax),fontsize=15,color='k')
    
cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
    
plt.savefig(directoryfigure + 'TrendsMap_%s_%s-%s_%s.png' % (sliceperiod,yearmin,yearmax,variq),dpi=300)
