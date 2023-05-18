"""
Plot trends in FMA from 1979 to 2022 for SND
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
directorydata = '/work/Zachary.Labe/Data/ERA5_MEDS/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'FMA'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True

### Read data
data = Dataset(directorydata + 'SND_1979-2021.nc')
lat1 = data.variables['lat'][:]
lon1 = data.variables['lon'][:]
data.set_auto_mask(False) # weird fill values after regridding
snd = data.variables['SND'][:].reshape(years.shape[0],12,lat1.shape[0],lon1.shape[0])
data.close()

if sliceperiod == 'April':
    datat = np.nanmean(snd[:,3:4,:,:],axis=1)
elif sliceperiod == 'FMA':
    datat = np.nanmean(snd[:,1:4,:,:],axis=1)

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2021
trend = UT.linearTrendR(datat,years,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-0.03,0.031,0.001)
barlim = np.round(np.arange(-0.03,0.031,0.03),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{Snow Depth Trend [m/decade]}'

ax = plt.subplot(111)

var = trend

m = Basemap(projection='ortho',lon_0=265,lat_0=60,resolution='l',area_thresh=10000)
m.drawcoastlines(color='darkgrey',linewidth=0.75)
# m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='darkgrey',linewidth=0.35)
    
# var, lons_cyclic = addcyclic(var, lon1)
# var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
# lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
# x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='k',
                  linewidth=0.7)
circle.set_clip_on(False)
m.drawlsmask(land_color=(0,0,0,0),ocean_color='k',lakes=False,zorder=11)

lon2,lat2 = np.meshgrid(lon1,lat1)
cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)

cs1.set_cmap(cmocean.cm.balance_r)

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

plt.title(r'\textbf{FEBRUARY -- APRIL}',fontsize=20,color='w')
    
cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='w',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('darkgrey')
plt.tight_layout()
    
plt.savefig(directoryfigure + 'PRESENT_TrendsMap_%s_%s-%s_SND.png' % (sliceperiod,yearmin,yearmax),dpi=300)
