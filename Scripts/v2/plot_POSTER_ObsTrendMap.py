"""
Plot trends in April from 1979 for T2M for AGU

Author    : Zachary M. Labe
Date      : 8 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthlyHighRes as ERA
import cmasher as cmr

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/AGU_Poster/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Parameters
variq = 'T2M'
sliceperiod = 'FMA'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,temp = ERA.read_ERA5_monthlyHighRes(variq,directorydata,sliceperiod,
                                                    years,sliceshape,addclimo,
                                                    slicenan)

### Read in 2022 data
data = Dataset('/work/Zachary.Labe/Data/ERA5_025x025/' + 'T2M_2022_Oct.nc')
nrtq = np.nanmean(data.variables['T2M'][:],axis=1) - 273.15
data.close()

### April
nrt = np.nanmean(nrtq[np.newaxis,1:4,:,:],axis=1)

### Combinate timeseries
yearsn = np.arange(1979,2022+1,1)
tempn = np.append(temp,nrt,axis=0)

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2022
trend = UT.linearTrendR(tempn,yearsn,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################
### Plot map of trends
def setcolor(x, color):
      for m in x:
          for t in x[m][1]:
              t.set_color(color)

limit = np.arange(-1,1.001,0.005)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{Trend [$^{\circ}$C/decade]}'

ax = plt.subplot(111)

var = trend

m = Basemap(projection='ortho',lon_0=255,lat_0=46,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
# m.drawstates(color='dimgrey',linewidth=0.3)
# m.drawcountries(color='dimgrey',linewidth=0.3)
    
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

parallels = np.arange(-90,91,30)
meridians = np.arange(-180,180,60)
par=m.drawparallels(parallels,labels=[True,True,True,True],linewidth=1,
                color='darkgrey',fontsize=4,zorder=6)
mer=m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=1,
                    fontsize=4,color='darkgrey',zorder=6)
setcolor(mer,'k')
setcolor(par,'k')

cs1 = m.contourf(x,y,var,limit,extend='both')

cs1.set_cmap(cmocean.cm.balance)

## Box 2
la1 = 43
la2 = 60
lo1 = 240
lo2 = 295
lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
latsslice = np.ones(len(lonsslice))*la2
m.plot(lonsslice, latsslice, color='gold', linewidth=5, latlon=True,zorder=7)
latsslice = np.ones(len(lonsslice))*la1
m.plot(lonsslice, latsslice, color='gold', linewidth=5, latlon=True,zorder=7)
m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=5,color='gold',zorder=7)
m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=5,color='gold',zorder=7)

cbar_ax1 = fig.add_axes([0.01,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=10,color='darkgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=10)
cbar1.outline.set_edgecolor('darkgrey')
plt.tight_layout()
    
plt.savefig(directoryfigure + 'POSTER_TrendsMap_%s_%s-%s_%s.png' % (sliceperiod,yearmin,yearmax,variq),dpi=600)
