"""
Plot differences in trends from observational datasets
Author    : Zachary M. Labe
Date      : 23 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import calc_Utilities as UT
import sys
import read_ERA5_monthly as ERA
import read_NCEP2 as NP
import read_BEST as B
import read_GISTEMP as G
import read_NOAAglobaltemp as NO
import scipy.stats as sts
import sys
import palettable.cubehelix as cm

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
datalabels = ['ERA5','NCEP2','GISTEMPv4','BEST','NOAA','OBS-MEAN']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
timelabels = ['1979-2020','1979-2020','1979-2020','1979-2020','1979-2020','1979-2020',
              '1979-2010','1979-2010','1979-2010','1979-2010','1979-2010','1979-2010']
variq = 'T2M'
sliceperiod = 'April'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
alpha = 0.1

### Read data
lat1,lon1,era5 = ERA.read_ERA5_monthly(variq,'/work/Zachary.Labe/Data/ERA5_19x25/',sliceperiod,years,sliceshape,addclimo,slicenan,False)
lat1,lon1,ncep = NP.read_NCEP2('/work/Zachary.Labe/Data/NCEP2/',sliceperiod,years,sliceshape,addclimo,slicenan,False)
lat1,lon1,giss = G.read_GISTEMP('/work/Zachary.Labe/Data/GISTEMP/',sliceperiod,years,sliceshape,addclimo,slicenan,False)
lat1b,lon1b,best = B.read_BEST('/work/Zachary.Labe/Data/BEST/',sliceperiod,years,sliceshape,addclimo,slicenan)
lat1,lon1,noaa = NO.read_NOAA('/work/Zachary.Labe/Data/NOAAGlobalTemp/',sliceperiod,years,sliceshape,addclimo,slicenan,False)

trende20 = UT.linearTrendR(np.asarray(era5),years,'surface',1979,2020)*10
trendn20 = UT.linearTrendR(np.asarray(ncep),years,'surface',1979,2020)*10
trendg20 = UT.linearTrendR(np.asarray(giss),years,'surface',1979,2020)*10
trendb20 = UT.linearTrendR(np.asarray(best),years,'surface',1979,2020)*10
trendno20 = UT.linearTrendR(np.asarray(noaa),years,'surface',1979,2020)*10

trende10 = UT.linearTrendR(np.asarray(era5),years,'surface',1979,2010)*10
trendn10 = UT.linearTrendR(np.asarray(ncep),years,'surface',1979,2010)*10
trendg10 = UT.linearTrendR(np.asarray(giss),years,'surface',1979,2010)*10
trendb10 = UT.linearTrendR(np.asarray(best),years,'surface',1979,2010)*10
trendno10 = UT.linearTrendR(np.asarray(noaa),years,'surface',1979,2010)*10

### Calculate statistical tests for 1979-2020
yearq20 = np.where((years == 2020))[0][0]+1
pvalse20 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsn20 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsg20 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsb20 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsno20 = np.empty((lat1.shape[0],lon1.shape[0]))
for i in range(lat1.shape[0]):
    for j in range(lon1.shape[0]):
        trenddefe20,he20,pvalse20[i,j],z = UT.mk_test(era5[:yearq20,i,j],alpha)
        trenddefn20,hn20,pvalsn20[i,j],z = UT.mk_test(ncep[:yearq20,i,j],alpha)
        trenddefg20,hg20,pvalsg20[i,j],z = UT.mk_test(giss[:yearq20,i,j],alpha)
        trenddefb20,hb20,pvalsb20[i,j],z = UT.mk_test(best[:yearq20,i,j],alpha)
        trenddefno20,hno20,pvalsno20[i,j],z = UT.mk_test(noaa[:yearq20,i,j],alpha)
        
pvalse20[np.where(pvalse20 == 1.)] = 0.
pvalse20[np.where(np.isnan(pvalse20))] = 1.
pvalse20[np.where(pvalse20 == 0.)] = np.nan

pvalsn20[np.where(pvalsn20 == 1.)] = 0.
pvalsn20[np.where(np.isnan(pvalsn20))] = 1.
pvalsn20[np.where(pvalsn20 == 0.)] = np.nan

pvalsg20[np.where(pvalsg20 == 1.)] = 0.
pvalsg20[np.where(np.isnan(pvalsg20))] = 1.
pvalsg20[np.where(pvalsg20 == 0.)] = np.nan

pvalsb20[np.where(pvalsb20 == 1.)] = 0.
pvalsb20[np.where(np.isnan(pvalsb20))] = 1.
pvalsb20[np.where(pvalsb20 == 0.)] = np.nan

pvalsno20[np.where(pvalsno20 == 1.)] = 0.
pvalsno20[np.where(np.isnan(pvalsno20))] = 1.
pvalsno20[np.where(pvalsno20 == 0.)] = np.nan
        
### Calculate statistical tests for 1979-2010
yearq10 = np.where((years == 2010))[0][0]+1
pvalse10 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsn10 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsg10 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsb10 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsno10 = np.empty((lat1.shape[0],lon1.shape[0]))
for i in range(lat1.shape[0]):
    for j in range(lon1.shape[0]):
        trenddefe10,he10,pvalse10[i,j],z = UT.mk_test(era5[:yearq10,i,j],alpha)
        trenddefn10,hn10,pvalsn10[i,j],z = UT.mk_test(ncep[:yearq10,i,j],alpha)
        trenddefg10,hg10,pvalsg10[i,j],z = UT.mk_test(giss[:yearq10,i,j],alpha)
        trenddefb10,hb10,pvalsb10[i,j],z = UT.mk_test(best[:yearq10,i,j],alpha)
        trenddefno10,hno10,pvalsno10[i,j],z = UT.mk_test(noaa[:yearq10,i,j],alpha)
        
pvalse10[np.where(pvalse10 == 1.)] = 0.
pvalse10[np.where(np.isnan(pvalse10))] = 1.
pvalse10[np.where(pvalse10 == 0.)] = np.nan

pvalsn10[np.where(pvalsn10 == 1.)] = 0.
pvalsn10[np.where(np.isnan(pvalsn10))] = 1.
pvalsn10[np.where(pvalsn10 == 0.)] = np.nan

pvalsg10[np.where(pvalsg10 == 1.)] = 0.
pvalsg10[np.where(np.isnan(pvalsg10))] = 1.
pvalsg10[np.where(pvalsg10 == 0.)] = np.nan

pvalsb10[np.where(pvalsb10 == 1.)] = 0.
pvalsb10[np.where(np.isnan(pvalsb10))] = 1.
pvalsb10[np.where(pvalsb10 == 0.)] = np.nan

pvalsno10[np.where(pvalsno10 == 1.)] = 0.
pvalsno10[np.where(np.isnan(pvalsno10))] = 1.
pvalsno10[np.where(pvalsno10 == 0.)] = np.nan

### Calculate observational means
mean20 = np.nanmean(np.array([trende20,trendn20,trendg20,trendb20,trendno20]),axis=0)
mean10 = np.nanmean(np.array([trende10,trendn10,trendg10,trendb10,trendno10]),axis=0)
emptypval20 = np.full(mean20.shape,np.nan)
emptypval10 = np.full(mean10.shape,np.nan)

### Create lists for plotting
plottrend20 = [trende20,trendn20,trendg20,trendb20,trendno20,mean20] 
plottrend10 = [trende10,trendn10,trendg10,trendb10,trendno10,mean10] 
plotpval20 = [pvalse20,pvalsn20,pvalsg20,pvalsb20,pvalsno20,emptypval20] 
plotpval10 = [pvalse10,pvalsn10,pvalsg10,pvalsb10,pvalsno10,emptypval10] 

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of LRP means training
limit = np.arange(-0.75,0.751,0.025)
barlim = np.arange(-0.75,0.76,0.25)
cmap = cmocean.cm.balance
label = r'\textbf{T2M [$^{\circ}$C/decade]}'

fig = plt.figure(figsize=(10,3))
for r in range(len(plottrend20)*2):
    if r < 6:
        var = plottrend20[r]
        pvar = plotpval20[r]
    else:
        var = plottrend10[r-6]
        pvar = plotpval10[r-6]        
    
    ax1 = plt.subplot(2,6,r+1)
    m = Basemap(projection='ortho',lon_0=265,lat_0=40,resolution='l',area_thresh=10000)
    circle = m.drawmapboundary(fill_color='k')
    circle.set_clip_on(False) 
    m.drawcoastlines(color='dimgrey',linewidth=0.35)
    
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    pvar,lons_cyclic = addcyclic(pvar, lon1)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
    x, y = m(lon2d, lat2d)
       
    circle = m.drawmapboundary(fill_color='darkgray',color='darkgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs = m.contourf(x,y,var,limit,extend='both')
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['.....'])
    cs.set_cmap(cmap) 
    
    m.drawlsmask(land_color=(0,0,0,0),ocean_color='dimgrey',lakes=False,zorder=11)
    
    if r < 6:
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
        latsslice = np.ones(len(lonsslice))*la2
        m.plot(lonsslice, latsslice, color='gold', linewidth=0.8, latlon=True,zorder=4)
        latsslice = np.ones(len(lonsslice))*la1
        m.plot(lonsslice, latsslice, color='gold', linewidth=0.8, latlon=True,zorder=4)
        m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=0.8,color='gold',zorder=4)
        m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=0.8,color='gold',zorder=4)\
    
    if r < 6:
        plt.title(r'\textbf{%s}' % datalabels[r],color='dimgrey',fontsize=15)
         
    ax1.annotate(r'\textbf{[%s]}' % letters[r],xy=(0,0),xytext=(0.86,0.97),
                  textcoords='axes fraction',color='k',fontsize=6,
                  rotation=330,ha='center',va='center')
    
    if any([r==0,r==6]):
        ax1.annotate(r'\textbf{%s}' % timelabels[r],xy=(0,0),xytext=(-0.07,0.5),
                      textcoords='axes fraction',color='k',fontsize=10,
                      rotation=90,ha='center',va='center')
    
###############################################################################
cbar_ax1 = fig.add_axes([0.352,0.11,0.3,0.03])                
cbar1 = fig.colorbar(cs,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=10,color='dimgrey',labelpad=0.7)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')

plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=-0.1,bottom=0.17,top=0.90)
plt.savefig(directoryfigure + 'differencesObservations_MapsOfTrends_combinedTimes_%s.png' % sliceperiod,dpi=300)
