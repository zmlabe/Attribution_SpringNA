"""
Plot trends in FMA from 1979 to 2019 for Z500
Author    : Zachary M. Labe
Date      : 6 June 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthlyBE as ERA

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ERA5/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'Z500'
sliceperiod = 'FMA'
yearsBE = np.arange(1950,2019+1,1)
years = np.arange(1979,2019+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,data = ERA.read_ERA5_monthlyBE(variq,directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan)

### Slice out 1950-1979
yearq = np.where((yearsBE >= 1979))[0]
fma = data[yearq,:,:]

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2019
trend = UT.linearTrendR(fma,years,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-15,15.01,0.5)
barlim = np.round(np.arange(-15,16,5),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{Z500 Trend m/decade]}'

ax = plt.subplot(111)

var = trend

m = Basemap(projection='ortho',lon_0=265,lat_0=90,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)
    
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,var,limit,extend='both')

cs1.set_cmap(cmocean.cm.balance)

## Box 2
la1 = 43
la2 = 60
lo1 = 240
lo2 = 295
lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
latsslice = np.ones(len(lonsslice))*la2
m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
latsslice = np.ones(len(lonsslice))*la1
m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1.5,color='b',zorder=4)
m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1.5,color='b',zorder=4)

plt.title(r'\textbf{FEBRUARY -- APRIL %s-%s}' % (yearmin,yearmax),fontsize=15,color='k')
    
cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
    
plt.savefig(directoryfigure + 'TrendsMap_FMA_%s-%s_Z500.png' % (yearmin,yearmax),dpi=300)
