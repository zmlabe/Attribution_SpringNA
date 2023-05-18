"""
Plot trends in spring from 1979 for SF200

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
variq = 'SF200'
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

### Convert for plotting
data = data/1e6

### Calculate anomalies
anoms = calc_anomalies(years,data)

### Pick years between 1979 and 2020
yearq = np.where((years >= 1979) & (years <= 2020))[0]
anomnew = anoms[yearq]
datanew = data[yearq]

### Composite time periods
period1q = np.where((years >= 1980) & (years <= 1994))[0]
SF2001 = np.nanmean(anomnew[period1q,:,:],axis=0)
SF2001c = np.nanmean(datanew[period1q,:,:],axis=0)

period2q = np.where((years >= 1985) & (years <= 1999))[0]
SF2002 = np.nanmean(anomnew[period2q,:,:],axis=0)
SF2002c = np.nanmean(datanew[period2q,:,:],axis=0)

period3q = np.where((years >= 1990) & (years <= 2004))[0]
SF2003 = np.nanmean(anomnew[period3q,:,:],axis=0)
SF2003c = np.nanmean(datanew[period3q,:,:],axis=0)

period4q = np.where((years >= 1995) & (years <= 2009))[0]
SF2004 = np.nanmean(anomnew[period4q,:,:],axis=0)
SF2004c = np.nanmean(datanew[period4q,:,:],axis=0)

period5q = np.where((years >= 2000) & (years <= 2014))[0]
SF2005 = np.nanmean(anomnew[period5q,:,:],axis=0)
SF2005c = np.nanmean(datanew[period5q,:,:],axis=0)

period6q = np.where((years >= 2005) & (years <= 2019))[0]
SF2006 = np.nanmean(anomnew[period6q,:,:],axis=0)
SF2006c = np.nanmean(datanew[period6q,:,:],axis=0)

### 1981-2010 climatology
yearbaselineq = np.where((years >= 1981) & (years <= 2010))[0]
meanSF200c = np.nanmean(datanew[yearbaselineq,:,:],axis=0)

SF200_all = [SF2001,SF2002,SF2003,SF2004,SF2005,SF2006]
# SF200_allc = [SF2001c,SF2002c,SF2003c,SF2004c,SF2005c,SF2006c]
SF200_allc = [meanSF200c,meanSF200c,meanSF200c,meanSF200c,meanSF200c,meanSF200c]

###############################################################################
###############################################################################
###############################################################################
### Large pacific view

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-5.01,5.011,0.1)
barlim = np.round(np.arange(-5,6,5),2)
cmap = cmocean.cm.balance
label = r'\textbf{SF200 Anomalies [10$^{6}$ m$^{2}$ s$^{-1}$] -- ERA5}'

fig = plt.figure(figsize=(4,10))
for i in range(len(SF200_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=100, llcrnrlat=0,
                urcrnrlon=270, urcrnrlat=65,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,SF200_all[i],limit,extend='both',latlon=True)
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
     
    plt.text(76,5,r'\textbf{%s}' % yearslist[i],rotation=90,color='k',fontsize=14)   
  
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
                
plt.savefig(directoryfigure + 'SF200_%s_AllTimePeriods_Obs.png' % sliceperiod,dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'SF850'
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

### Convert for plotting
data = data/1e6

### Calculate anomalies
anoms = calc_anomalies(years,data)

### Pick years between 1979 and 2020
yearq = np.where((years >= 1979) & (years <= 2020))[0]
anomnew = anoms[yearq]
datanew = data[yearq]

### Composite time periods
period1q = np.where((years >= 1980) & (years <= 1994))[0]
SF8501 = np.nanmean(anomnew[period1q,:,:],axis=0)
SF8501c = np.nanmean(datanew[period1q,:,:],axis=0)

period2q = np.where((years >= 1985) & (years <= 1999))[0]
SF8502 = np.nanmean(anomnew[period2q,:,:],axis=0)
SF8502c = np.nanmean(datanew[period2q,:,:],axis=0)

period3q = np.where((years >= 1990) & (years <= 2004))[0]
SF8503 = np.nanmean(anomnew[period3q,:,:],axis=0)
SF8503c = np.nanmean(datanew[period3q,:,:],axis=0)

period4q = np.where((years >= 1995) & (years <= 2009))[0]
SF8504 = np.nanmean(anomnew[period4q,:,:],axis=0)
SF8504c = np.nanmean(datanew[period4q,:,:],axis=0)

period5q = np.where((years >= 2000) & (years <= 2014))[0]
SF8505 = np.nanmean(anomnew[period5q,:,:],axis=0)
SF8505c = np.nanmean(datanew[period5q,:,:],axis=0)

period6q = np.where((years >= 2005) & (years <= 2019))[0]
SF8506 = np.nanmean(anomnew[period6q,:,:],axis=0)
SF8506c = np.nanmean(datanew[period6q,:,:],axis=0)

### 1981-2010 climatology
yearbaselineq = np.where((years >= 1981) & (years <= 2010))[0]
meanSF850c = np.nanmean(datanew[yearbaselineq,:,:],axis=0)

SF850_all = [SF8501,SF8502,SF8503,SF8504,SF8505,SF8506]
# SF850_allc = [SF8501c,SF8502c,SF8503c,SF8504c,SF8505c,SF8506c]
SF850_allc = [meanSF850c,meanSF850c,meanSF850c,meanSF850c,meanSF850c,meanSF850c]

###############################################################################
###############################################################################
###############################################################################
### Large pacific view

yearslist = ['1980-1994','1985-1999','1990-8504','1995-8509','2000-2014','2005-2019']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-5.01,5.011,0.1)
barlim = np.round(np.arange(-5,6,5),2)
cmap = cmocean.cm.balance
label = r'\textbf{SF850 Anomalies [10$^{6}$ m$^{2}$ s$^{-1}$] -- ERA5}'

fig = plt.figure(figsize=(4,10))
for i in range(len(SF850_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=100, llcrnrlat=0,
                urcrnrlon=270, urcrnrlat=65,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,SF850_all[i],limit,extend='both',latlon=True)
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
     
    plt.text(76,5,r'\textbf{%s}' % yearslist[i],rotation=90,color='k',fontsize=14)   
  
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
                
plt.savefig(directoryfigure + 'SF850_%s_AllTimePeriods_Obs.png' % sliceperiod,dpi=300)
