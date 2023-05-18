"""
Plot trends in spring from 1979 for SST

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
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/Presentations/' 
directorydata = '/work/Zachary.Labe/Data/' 

def calc_anomalies(years,data):
    """ 
    Calculate anomalies
    """
    
    ### Baseline - 1981-2010
    if data.ndim == 3:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[yearqold,:,:],axis=0)
        anoms = data - climold
    elif data.ndim == 4:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:]
    
    return anoms

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'SST'
sliceperiod = 'FMA'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,lev1,data = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,'surface')

lon2,lat2 = np.meshgrid(lon1,lat1)

### Calculate anomalies
anoms = calc_anomalies(years,data)

### Pick years between 1979 and 2020
yearq = np.where((years >= 1979) & (years <= 2020))[0]
anomnew = anoms[yearq]

### Composite time periods
period1q = np.where((years >= 1980) & (years <= 1994))[0]
sst1 = np.nanmean(anomnew[period1q,:,:],axis=0)

period2q = np.where((years >= 1985) & (years <= 1999))[0]
sst2 = np.nanmean(anomnew[period2q,:,:],axis=0)

period3q = np.where((years >= 1990) & (years <= 2004))[0]
sst3 = np.nanmean(anomnew[period3q,:,:],axis=0)

period4q = np.where((years >= 1995) & (years <= 2009))[0]
sst4 = np.nanmean(anomnew[period4q,:,:],axis=0)

period5q = np.where((years >= 2000) & (years <= 2014))[0]
sst5 = np.nanmean(anomnew[period5q,:,:],axis=0)

period6q = np.where((years >= 2005) & (years <= 2019))[0]
sst6 = np.nanmean(anomnew[period6q,:,:],axis=0)

sst_all = [sst1,sst2,sst3,sst4,sst5,sst6]

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-1.01,1.011,0.001)
barlim = np.round(np.arange(-1,2,1),2)
cmap = cmocean.cm.balance
label = r'\textbf{SST Anomalies [$^{\circ}$C] -- ERA5}'

fig = plt.figure(figsize=(5,10))
for i in range(len(sst_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=120, llcrnrlat=20,
                urcrnrlon=240, urcrnrlat=60,resolution='l',area_thresh=1000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,sst_all[i],limit,extend='both',latlon=True)
    cs1.set_cmap(cmap) 
    
    if i == 0:
        parallels = np.arange(-90,91,20)
        meridians = np.arange(-180,180,20)
        par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                        color='w',fontsize=4,zorder=40)
        mer=m.drawmeridians(meridians,labels=[False,False,True,False],linewidth=0.5,
                            fontsize=4,color='w',zorder=40)
    elif any([i==1,i==2,i==3,i==4]):
        parallels = np.arange(-90,91,20)
        meridians = np.arange(-180,180,20)
        par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                        color='w',fontsize=4,zorder=40)
        mer=m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.5,
                            fontsize=4,color='w',zorder=40)
    elif i == 5:
        parallels = np.arange(-90,91,20)
        meridians = np.arange(-180,180,20)
        par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                        color='w',fontsize=4,zorder=40)
        mer=m.drawmeridians(meridians,labels=[False,False,False,True],linewidth=0.5,
                            fontsize=4,color='w',zorder=40)
     
    plt.text(104,26,r'\textbf{%s}' % yearslist[i],rotation=90,color='k',fontsize=14)   
  
###############################################################################
cbar_ax1 = fig.add_axes([0.37,0.06,0.3,0.02])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=10,color='dimgrey',labelpad=1.8)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')

# plt.tight_layout()
plt.subplots_adjust(hspace=-0.01)
                
plt.savefig(directoryfigure + 'PRESENT_SST_%s_AllTimePeriods_Obs.png' % sliceperiod,dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Large pacific view

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-1.01,1.011,0.001)
barlim = np.round(np.arange(-1,2,1),2)
cmap = cmocean.cm.balance
label = r'\textbf{SST Anomalies [$^{\circ}$C] -- ERA5}'

fig = plt.figure(figsize=(3.5,10))
for i in range(len(sst_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=100, llcrnrlat=-20,
                urcrnrlon=270, urcrnrlat=65,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,sst_all[i],limit,extend='both',latlon=True)
    cs1.set_cmap(cmap) 
    
    # if i == 0:
    #     parallels = np.arange(-90,91,20)
    #     meridians = np.arange(-180,180,20)
    #     par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
    #                     color='w',fontsize=4,zorder=40)
    #     mer=m.drawmeridians(meridians,labels=[False,False,True,False],linewidth=0.5,
    #                         fontsize=4,color='w',zorder=40)
    # elif any([i==1,i==2,i==3,i==4]):
    #     parallels = np.arange(-90,91,20)
    #     meridians = np.arange(-180,180,20)
    #     par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
    #                     color='w',fontsize=4,zorder=40)
    #     mer=m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.5,
    #                         fontsize=4,color='w',zorder=40)
    # elif i == 5:
    #     parallels = np.arange(-90,91,20)
    #     meridians = np.arange(-180,180,20)
    #     par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
    #                     color='w',fontsize=4,zorder=40)
    #     mer=m.drawmeridians(meridians,labels=[False,False,False,True],linewidth=0.5,
    #                         fontsize=4,color='w',zorder=40)
     
    plt.text(70,-7,r'\textbf{%s}' % yearslist[i],rotation=90,color='w',fontsize=14)   
  
###############################################################################
cbar_ax1 = fig.add_axes([0.37,0.06,0.3,0.02])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=10,color='darkgrey',labelpad=1.8)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('darkgrey')

# plt.tight_layout()
plt.subplots_adjust(hspace=-0.01)
                
plt.savefig(directoryfigure + 'PRESENT_SST_%s_AllTimePeriods_Obs.png' % sliceperiod,dpi=600)
