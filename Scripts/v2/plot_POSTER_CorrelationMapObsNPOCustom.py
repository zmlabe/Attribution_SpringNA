"""
Create a custom index for T2M_BoxNA

Author    : Zachary M. Labe
Date      : 30 November 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
import calc_DetrendData as DT
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import cmocean
import cmasher as cmr

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/AGU_Poster/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/work/Zachary.Labe/Data/ClimateIndices/'
directorydataoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/Obs/'

### Select Parameters
variqM = 'Z500'
sliceperiod = 'April'
yearmin = 1979
yearmax = 2020
years = np.arange(yearmin,yearmax+1,1)

### Read data
corrmaps = np.load(directorydataoutput + 'ObsCorrelationPatterns_T2M_BoxNApattern_%s-%s_%s-%s.npz' % (variqM,sliceperiod,'T2M_BoxNA',sliceperiod))
lat = corrmaps['lat']
lon = corrmaps['lon']
lon2,lat2 = np.meshgrid(lon,lat)

pvalue = corrmaps['pval']
rvalue = corrmaps['r']

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-1,1.01,0.1)
barlim = np.round(np.arange(-1,1.1,0.5),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{Correlation Coefficient}'

ax = plt.subplot(111)

var = rvalue
pvar = pvalue

pvar[np.isnan(pvar)] = 0.
var = var*pvar
var[np.where(var == 0.)] = np.nan

m = Basemap(projection='ortho',lon_0=230,lat_0=55,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1.5)
# m.drawstates(color='dimgrey',linewidth=0.3)
m.drawcountries(color='dimgrey',linewidth=0.6)
    
pvar,lons_cyclic = addcyclic(pvar, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
var, lons_cyclic = addcyclic(var, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
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

cs1 = m.contourf(x,y,var,limit,extend='both')

cs1.set_cmap(cmr.fusion_r)

## Box 3 - Custom (North)
la1 = 55
la2 = 70
lo1 = 175
lo2 = 205
lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
latsslice = np.ones(len(lonsslice))*la2
m.plot(lonsslice, latsslice, color='k', linewidth=4, latlon=True,zorder=4)
latsslice = np.ones(len(lonsslice))*la1
m.plot(lonsslice, latsslice, color='k', linewidth=4, latlon=True,zorder=4)
m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=4,color='k',zorder=4)
m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=4,color='k',zorder=4)

## Box 4 - Custom (South)
la1 = 22
la2 = 35
lo1 = 183
lo2 = 213
lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
latsslice = np.ones(len(lonsslice))*la2
m.plot(lonsslice, latsslice, color='k', linewidth=4, latlon=True,zorder=4)
latsslice = np.ones(len(lonsslice))*la1
m.plot(lonsslice, latsslice, color='k', linewidth=4, latlon=True,zorder=4)
m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=4,color='k',zorder=4)
m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=4,color='k',zorder=4)
   
cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=8,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()

plt.savefig(directoryfigure + 'POSTER_ObsCorrelationPatterns_T2M_BoxNApattern_%s-%s_%s-%s.png' % (variqM,sliceperiod,'T2M_BoxNA',sliceperiod),dpi=600)
