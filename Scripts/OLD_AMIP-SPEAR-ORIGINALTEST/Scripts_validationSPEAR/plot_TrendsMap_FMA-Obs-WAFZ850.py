"""
Plot trends in FMA from 1979 to 2022 for WAFz850
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
import read_ERA5_monthlyHighRes as ERA

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ERA5_1x1/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
sliceperiod = 'FMA'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True

### Read data
data = Dataset(directorydata + 'WAFZ850_1979-2021.nc')
lat1 = data.variables['latitude'][:]
lon1 = data.variables['longitude'][:]
snd = data.variables['WAFZ850'][:].reshape(years.shape[0],12,lat1.shape[0],lon1.shape[0])
data.close()

if sliceperiod == 'April':
    datat = np.asarray(np.nanmean(snd[:,3:4,:,:],axis=1))
    climo = np.nanmean(datat,axis=0)
elif sliceperiod == 'FMA':
    datat = np.asarray(np.nanmean(snd[:,1:4,:,:],axis=1))
    climo = np.nanmean(datat,axis=0)

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2019
trend = UT.linearTrendR(datat,years,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-0.02,0.021,0.001)
barlim = np.round(np.arange(-0.02,0.04,0.02),4)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{WAFZ850 Trend ["m$^{2}$/s$^{2}$/decade]}'

ax = plt.subplot(111)

var = trend
cl = climo

m = Basemap(projection='npstere',lon_0=270,boundinglat=45,resolution='l',area_thresh=10000,
            round=True)
m.drawcoastlines(color='k',linewidth=1)
    
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
cl, lons_cyclic = addcyclic(cl, lon1)
cl, lons_cyclic = shiftgrid(180., cl, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,var,limit,extend='both')
# cs2 = m.contour(x,y,cl,np.arange(-0.3,0.31,0.05),extend='both',
#                 linewidths=2,colors='dimgrey')

cs1.set_cmap(cmocean.cm.balance_r)

## Box 2
# la1 = 43
# la2 = 60
# lo1 = 240
# lo2 = 295
# lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
# latsslice = np.ones(len(lonsslice))*la2
# m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
# latsslice = np.ones(len(lonsslice))*la1
# m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
# m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1.5,color='b',zorder=4)
# m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1.5,color='b',zorder=4)

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
    
plt.savefig(directoryfigure + 'TrendsMap_%s_%s-%s_WAFZ850.png' % (sliceperiod,yearmin,yearmax),dpi=300)
