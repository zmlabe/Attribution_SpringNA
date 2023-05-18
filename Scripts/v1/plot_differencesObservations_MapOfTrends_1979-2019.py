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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/ObsTrends/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
datalabels = ['ERA5','NCEP2','GISTEMPv4','BEST','NOAA','OBS-MEAN']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
timelabels = ['1979-2019','1979-2019','1979-2019','1979-2019','1979-2019','1979-2019',
              '1979-2010','1979-2010','1979-2010','1979-2010','1979-2010','1979-2010']
variq = 'T2M'
sliceperiod = 'MAM'
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

trende19 = UT.linearTrendR(np.asarray(era5),years,'surface',1979,2019)*10
trendn19 = UT.linearTrendR(np.asarray(ncep),years,'surface',1979,2019)*10
trendg19 = UT.linearTrendR(np.asarray(giss),years,'surface',1979,2019)*10
trendb19 = UT.linearTrendR(np.asarray(best),years,'surface',1979,2019)*10
trendno19 = UT.linearTrendR(np.asarray(noaa),years,'surface',1979,2019)*10

trende10 = UT.linearTrendR(np.asarray(era5),years,'surface',1979,2010)*10
trendn10 = UT.linearTrendR(np.asarray(ncep),years,'surface',1979,2010)*10
trendg10 = UT.linearTrendR(np.asarray(giss),years,'surface',1979,2010)*10
trendb10 = UT.linearTrendR(np.asarray(best),years,'surface',1979,2010)*10
trendno10 = UT.linearTrendR(np.asarray(noaa),years,'surface',1979,2010)*10

### Calculate statistical tests for 1979-2019
yearq19 = np.where((years == 2019))[0][0]+1
pvalse19 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsn19 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsg19 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsb19 = np.empty((lat1.shape[0],lon1.shape[0]))
pvalsno19 = np.empty((lat1.shape[0],lon1.shape[0]))
for i in range(lat1.shape[0]):
    for j in range(lon1.shape[0]):
        trenddefe19,he19,pvalse19[i,j],z = UT.mk_test(era5[:yearq19,i,j],alpha)
        trenddefn19,hn19,pvalsn19[i,j],z = UT.mk_test(ncep[:yearq19,i,j],alpha)
        trenddefg19,hg19,pvalsg19[i,j],z = UT.mk_test(giss[:yearq19,i,j],alpha)
        trenddefb19,hb19,pvalsb19[i,j],z = UT.mk_test(best[:yearq19,i,j],alpha)
        trenddefno19,hno19,pvalsno19[i,j],z = UT.mk_test(noaa[:yearq19,i,j],alpha)
        
pvalse19[np.where(pvalse19 == 1.)] = 0.
pvalse19[np.where(np.isnan(pvalse19))] = 1.
pvalse19[np.where(pvalse19 == 0.)] = np.nan

pvalsn19[np.where(pvalsn19 == 1.)] = 0.
pvalsn19[np.where(np.isnan(pvalsn19))] = 1.
pvalsn19[np.where(pvalsn19 == 0.)] = np.nan

pvalsg19[np.where(pvalsg19 == 1.)] = 0.
pvalsg19[np.where(np.isnan(pvalsg19))] = 1.
pvalsg19[np.where(pvalsg19 == 0.)] = np.nan

pvalsb19[np.where(pvalsb19 == 1.)] = 0.
pvalsb19[np.where(np.isnan(pvalsb19))] = 1.
pvalsb19[np.where(pvalsb19 == 0.)] = np.nan

pvalsno19[np.where(pvalsno19 == 1.)] = 0.
pvalsno19[np.where(np.isnan(pvalsno19))] = 1.
pvalsno19[np.where(pvalsno19 == 0.)] = np.nan
        
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
mean19 = np.nanmean(np.array([trende19,trendn19,trendg19,trendb19,trendno19]),axis=0)
mean10 = np.nanmean(np.array([trende10,trendn10,trendg10,trendb10,trendno10]),axis=0)
emptypval19 = np.full(mean19.shape,np.nan)
emptypval10 = np.full(mean10.shape,np.nan)

### Create lists for plotting
plottrend19 = [trende19,trendn19,trendg19,trendb19,trendno19,mean19] 
plottrend10 = [trende10,trendn10,trendg10,trendb10,trendno10,mean10] 
plotpval19 = [pvalse19,pvalsn19,pvalsg19,pvalsb19,pvalsno19,emptypval19] 
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
for r in range(len(plottrend19)*2):
    if r < 6:
        var = plottrend19[r]
        pvar = plotpval19[r]
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
plt.savefig(directoryfigure + 'differencesObservations_MapsOfTrends_combinedTimes13_%s.png' % sliceperiod,dpi=300)
