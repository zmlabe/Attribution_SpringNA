"""
Plot differences of climate baselines 1981-2010 and 1991-2020 by month
Author    : Zachary M. Labe
Date      : 16 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthly as ERA

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/Users/zlabe/Documents/Research/Attribution_SpringNA/Figures/' 
directorydata = '/Users/zlabe/Data/ERA5/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'none'
years = np.arange(1979,2022+1,1)
sliceshape = 4
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,data = ERA.read_ERA5_monthly(variq,directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,newon)

### Baseline - 1981-2010
yearqold = np.where((years >= 1981) & (years <= 2010))[0]
climold = np.nanmean(data[yearqold,:,:,:],axis=0)

### Baseline - 1991-2020
yearqold = np.where((years >= 1991) & (years <= 2020))[0]
climnew = np.nanmean(data[yearqold,:,:,:],axis=0)

### Difference
diff = climnew - climold

### February to April mean
diffmean = np.nanmean(diff[1:4,:,:],axis=0)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-1,1.01,0.02)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{1991-2020 minus 1981-2010 [$^{\circ}$]}'

ax = plt.subplot(111)

var = diffmean

m = Basemap(projection='ortho',lon_0=265,lat_0=40,resolution='l',area_thresh=10000)
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
la2 = 63
lo1 = 238
lo2 = 270
lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
latsslice = np.ones(len(lonsslice))*la2
m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
latsslice = np.ones(len(lonsslice))*la1
m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1.5,color='b',zorder=4)
m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1.5,color='b',zorder=4)

plt.title(r'\textbf{FEBRUARY -- APRIL}',fontsize=15,color='k')
    
cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
    
plt.savefig(directoryfigure + 'BaselineMapDifferences_FMA.png',dpi=300)
