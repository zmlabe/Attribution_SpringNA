"""
Plot trends in spring from 1979 for stationary waves

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
variq = 'Z300xwave1'
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
datanew = data[yearq]

### Composite time periods
period1q = np.where((years >= 1980) & (years <= 1994))[0]
z3001 = np.nanmean(anomnew[period1q,:,:],axis=0)
z3001c = np.nanmean(datanew[period1q,:,:],axis=0)

period2q = np.where((years >= 1985) & (years <= 1999))[0]
z3002 = np.nanmean(anomnew[period2q,:,:],axis=0)
z3002c = np.nanmean(datanew[period2q,:,:],axis=0)

period3q = np.where((years >= 1990) & (years <= 2004))[0]
z3003 = np.nanmean(anomnew[period3q,:,:],axis=0)
z3003c = np.nanmean(datanew[period3q,:,:],axis=0)

period4q = np.where((years >= 1995) & (years <= 2009))[0]
z3004 = np.nanmean(anomnew[period4q,:,:],axis=0)
z3004c = np.nanmean(datanew[period4q,:,:],axis=0)

period5q = np.where((years >= 2000) & (years <= 2014))[0]
z3005 = np.nanmean(anomnew[period5q,:,:],axis=0)
z3005c = np.nanmean(datanew[period5q,:,:],axis=0)

period6q = np.where((years >= 2005) & (years <= 2019))[0]
z3006 = np.nanmean(anomnew[period6q,:,:],axis=0)
z3006c = np.nanmean(datanew[period6q,:,:],axis=0)

### 1981-2010 climatology
yearbaselineq = np.where((years >= 1981) & (years <= 2010))[0]
meanz300c = np.nanmean(datanew[yearbaselineq,:,:],axis=0)

z300_all = [z3001,z3002,z3003,z3004,z3005,z3006]
# z300_allc = [z3001c,z3002c,z3003c,z3004c,z3005c,z3006c]
z300_allc = [meanz300c,meanz300c,meanz300c,meanz300c,meanz300c,meanz300c]

###############################################################################
###############################################################################
###############################################################################
### Large pacific view

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-20.01,20.011,0.01)
barlim = np.round(np.arange(-20,21,10),2)
cmap = cmocean.cm.balance
label = r'\textbf{z300x Anomalies [m] -- ERA5}'

fig = plt.figure(figsize=(4,10))
for i in range(len(z300_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=100, llcrnrlat=20,
                urcrnrlon=300, urcrnrlat=90,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,z300_all[i],limit,extend='both',latlon=True)
    cs2 = m.contour(lon2,lat2,z300_allc[i],np.arange(-200,201,20),extend='both',latlon=True,
                    linewidths=1,colors='k')
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
     
    plt.text(70,25,r'\textbf{%s}' % yearslist[i],rotation=90,color='k',fontsize=14)   
  
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
plt.subplots_adjust(hspace=-0.03)
                
plt.savefig(directoryfigure + 'z300x_%s_AllTimePeriods_Obs.png' % sliceperiod,dpi=300)

###############################################################################
###############################################################################
###############################################################################

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'Z500xwave1'
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
datanew = data[yearq]

### Composite time periods
period1q = np.where((years >= 1980) & (years <= 1994))[0]
z5001 = np.nanmean(anomnew[period1q,:,:],axis=0)
z5001c = np.nanmean(datanew[period1q,:,:],axis=0)

period2q = np.where((years >= 1985) & (years <= 1999))[0]
z5002 = np.nanmean(anomnew[period2q,:,:],axis=0)
z5002c = np.nanmean(datanew[period2q,:,:],axis=0)

period3q = np.where((years >= 1990) & (years <= 2004))[0]
z5003 = np.nanmean(anomnew[period3q,:,:],axis=0)
z5003c = np.nanmean(datanew[period3q,:,:],axis=0)

period4q = np.where((years >= 1995) & (years <= 2009))[0]
z5004 = np.nanmean(anomnew[period4q,:,:],axis=0)
z5004c = np.nanmean(datanew[period4q,:,:],axis=0)

period5q = np.where((years >= 2000) & (years <= 2014))[0]
z5005 = np.nanmean(anomnew[period5q,:,:],axis=0)
z5005c = np.nanmean(datanew[period5q,:,:],axis=0)

period6q = np.where((years >= 2005) & (years <= 2019))[0]
z5006 = np.nanmean(anomnew[period6q,:,:],axis=0)
z5006c = np.nanmean(datanew[period6q,:,:],axis=0)

### 1981-2010 climatology
yearbaselineq = np.where((years >= 1981) & (years <= 2010))[0]
meanz500c = np.nanmean(datanew[yearbaselineq,:,:],axis=0)

z500_all = [z5001,z5002,z5003,z5004,z5005,z5006]
# z500_allc = [z5001c,z5002c,z5003c,z5004c,z5005c,z5006c]
z500_allc = [meanz500c,meanz500c,meanz500c,meanz500c,meanz500c,meanz500c]

###############################################################################
###############################################################################
###############################################################################
### Large pacific view

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-20.01,20.011,0.01)
barlim = np.round(np.arange(-20,21,10),2)
cmap = cmocean.cm.balance
label = r'\textbf{z500x Anomalies [m] -- ERA5}'

fig = plt.figure(figsize=(4,10))
for i in range(len(z500_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=100, llcrnrlat=20,
                urcrnrlon=300, urcrnrlat=90,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,z500_all[i],limit,extend='both',latlon=True)
    cs2 = m.contour(lon2,lat2,z500_allc[i],np.arange(-200,201,20),extend='both',latlon=True,
                    linewidths=1,colors='k')
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
     
    plt.text(70,25,r'\textbf{%s}' % yearslist[i],rotation=90,color='k',fontsize=14)   
  
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
plt.subplots_adjust(hspace=-0.03)
                
plt.savefig(directoryfigure + 'z500x_%s_AllTimePeriods_Obs.png' % sliceperiod,dpi=300)
