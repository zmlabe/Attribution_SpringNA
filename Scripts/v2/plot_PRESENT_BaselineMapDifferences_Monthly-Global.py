"""
Plot differences of climate baselines 1981-2010 and 1991-2020 by month for 
global maps

Author    : Zachary M. Labe
Date      : 24 May 2022
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
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/Presentations/' 
directorydata = '/work/Zachary.Labe/Data/'

def setcolor(x, color):
      for m in x:
          for t in x[m][1]:
              t.set_color(color) 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'none'
years = np.arange(1979,2022+1,1)
sliceshape = 4
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,data = ERA.read_ERA5_monthlyHighRes(variq,directorydata,sliceperiod,
                                years,sliceshape,addclimo,slicenan)

### Baseline - 1981-2010
yearqold = np.where((years >= 1981) & (years <= 2010))[0]
climold = np.nanmean(data[yearqold,:,:,:],axis=0)

### Baseline - 1991-2020
yearqold = np.where((years >= 1991) & (years <= 2020))[0]
climnew = np.nanmean(data[yearqold,:,:,:],axis=0)

### Difference
diff = climnew - climold

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-1,1.01,0.05)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(6,6))
for mo in range(len(monthq)):                                                                                                                         
    label = r'\textbf{1991-2020 minus 1981-2010 [$^{\circ}$C]}'
    
    ax = plt.subplot(4,3,mo+1)
    
    var = diff[mo]

    m = Basemap(projection='robin',lon_0=0,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=0.3)
        
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
    x, y = m(lon2d, lat2d)
    
    circle = m.drawmapboundary(fill_color='k',color='k',linewidth=0.7)
    circle.set_clip_on(False)
    
    m.drawlsmask(land_color=(0,0,0,0),ocean_color='dimgrey',lakes=False,zorder=11)
    
    cs1 = m.contourf(x,y,var,limit,extend='both')
    
    cs1.set_cmap(cmocean.cm.balance)
    
    ax.annotate(r'\textbf{%s}' % (monthq[mo]),xy=(0,0),xytext=(0.03,0.90),
              textcoords='axes fraction',color='w',fontsize=9,
              rotation=40,ha='center',va='center')
    
cbar_ax1 = fig.add_axes([0.415,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=14,color='w',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('darkgrey')
plt.tight_layout()
plt.subplots_adjust(bottom=0.07,hspace=0.0,wspace=0.1)
    
plt.savefig(directoryfigure + 'PRESENT_BaselineMapDifferences_Monthly_Global.png',dpi=600)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-1,1.01,0.01)
barlim = np.round(np.arange(-1,2,1),2)

for mo in range(diff.shape[0]):                                                                                                                         
    label = r'\textbf{1991-2020 minus 1981-2010 [$^{\circ}$C]}'
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    var = diff[mo]

    m = Basemap(projection='robin',lon_0=0,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=0.5)
        
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
    x, y = m(lon2d, lat2d)
    
    circle = m.drawmapboundary(fill_color='k',color='k',linewidth=0.7)
    circle.set_clip_on(False)
    
    m.drawlsmask(land_color=(0,0,0,0),ocean_color='dimgrey',lakes=False,zorder=11)
    
    cs1 = m.contourf(x,y,var,limit,extend='both')
    
    cs1.set_cmap(cmocean.cm.balance)
    
    parallels = np.arange(-90,91,30)
    meridians = np.arange(-180,180,60)
    par=m.drawparallels(parallels,labels=[True,True,True,True],linewidth=0.6,
                    color='w',fontsize=0,zorder=12,clip_on=False)
    mer=m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=0.6,
                        fontsize=0,color='w',zorder=12,clip_on=False)
    setcolor(mer,'k')
    setcolor(par,'k')
    
    ax.annotate(r'\textbf{%s}' % (monthq[mo]),xy=(0,0),xytext=(0.06,0.90),
              textcoords='axes fraction',color='w',fontsize=20,
              rotation=40,ha='center',va='center')
    
    cbar_ax1 = fig.add_axes([0.43,0.08,0.15,0.025])                
    cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                        extend='both',extendfrac=0.07,drawedges=False)
    cbar1.set_label(label,fontsize=10,color='w',labelpad=1.4)  
    cbar1.set_ticks(barlim)
    cbar1.set_ticklabels(list(map(str,barlim)))
    cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
    cbar1.outline.set_edgecolor('darkgrey')
    plt.tight_layout()
    # plt.subplots_adjust(bottom=0.07,hspace=0.0,wspace=0.1)
        
    plt.savefig(directoryfigure + 'PRESENT_BaselineMapDifferences_%s_Global.png' % (monthq[mo]),dpi=600)
