"""
Plot maps of trends in slp and t2m for obs

Author    : Zachary M. Labe
Date      : 29 July 2022
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
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/Presentations/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
sliceperiod = 'April'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,lev1,t2m = ERA.read_ERA5_monthly1x1('T2M',directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,'surface')
lat1,lon1,lev1,slp = ERA.read_ERA5_monthly1x1('SLP',directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,'surface')

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2020
trendt2m = UT.linearTrendR(t2m,years,'surface',yearmin,yearmax)*10
trendslp = UT.linearTrendR(slp,years,'surface',yearmin,yearmax)*10

### Calculate climo
yearq = np.where((years >= 1981) & (years <= 2010))[0]
climo = np.nanmean(slp[yearq,:,:],axis=0)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different observational trends
limit = np.arange(-0.75,0.751,0.01)
barlim = np.arange(-0.75,0.76,0.75)

fig = plt.figure(figsize=(8,4))
label = r'\textbf{T2M Trend [$^{\circ}$C/decade]}'

ax = plt.subplot(1,2,1)

var = trendt2m

m = Basemap(projection='ortho',lon_0=230,lat_0=40,resolution='l',area_thresh=10000)
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
m.plot(lonsslice, latsslice, color='gold', linewidth=3, latlon=True,zorder=4)
latsslice = np.ones(len(lonsslice))*la1
m.plot(lonsslice, latsslice, color='gold', linewidth=3, latlon=True,zorder=4)
m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=3,color='gold',zorder=4)
m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=3,color='gold',zorder=4)

ax.annotate(r'\textbf{[a]}',xy=(0,0),xytext=(0.86,0.92),
              textcoords='axes fraction',color='k',fontsize=10,
              rotation=0,ha='center',va='center')
    
cbar_ax1 = fig.add_axes([0.02,0.08,0.1,0.02])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=6,color='darkgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=6)
cbar1.outline.set_edgecolor('darkgrey')

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different observational trends
limit = np.arange(-2,2.1,0.1)
barlim = np.arange(-2,3,2)
label = r'\textbf{SLP Trend [hPa/decade]}'

ax = plt.subplot(1,2,2)

var = trendslp

m = Basemap(projection='ortho',lon_0=230,lat_0=40,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)
    
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
climo, lons_cyclic = addcyclic(climo, lon1)
climo, lons_cyclic = shiftgrid(180., climo, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,var,limit,extend='both')
cs2 = m.contour(x,y,climo,np.arange(960,1050,3),extend='both',
                linewidths=1.2,colors='k')

cs1.set_cmap(cmocean.cm.balance)

ax.annotate(r'\textbf{[b]}',xy=(0,0),xytext=(0.86,0.92),
              textcoords='axes fraction',color='k',fontsize=10,
              rotation=0,ha='center',va='center')
    
cbar_ax1 = fig.add_axes([0.88,0.08,0.1,0.02])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=6,color='darkgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=6)
cbar1.outline.set_edgecolor('darkgrey')

plt.tight_layout()
    
plt.savefig(directoryfigure + 'PRESENT_TrendsMap_Obs_%s_%s-%s_SLPT2M.png' % (sliceperiod,yearmin,yearmax),dpi=300)
