"""
Compare trends in the data for the AMIPS from 1979 - 2019 for only SPEAR
along with its coupled run
 
Author    : Zachary M. Labe
Date      : 30 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools
import read_ERA5_monthlyLOWS_hindcastSPEAR_65lev as ER
import read_SPEAR_LOW_Hindcast_33lev as SPlow
import read_SPEAR_LOW_Hindcast_65lev as SPhigh
import read_SPEAR_MED as SPM
import read_FACTS_AMIPS as FA

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/SPEAR_Low_65lev/' 
directorydata = '/work/Zachary.Labe/Data/SPEAR/Hindcasts/SPEAR_LOW/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
scenario = ['amip_1880s_rf','amip_obs_rf','spear']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
experi1 = ['ERA5','SPEAR-LOW-65Lev','SPEAR-LOW-33Lev','SPEAR-MED-AMIP','SPEAR-MED']
model = [model_1880s_rf,model_obs_rf,model_spear]
modelunravel = list(itertools.chain(*model))
variq = 'T2M'
slicemonth = 'FMA'
monthhind = 'Feb'
slicenan = 'nan'
yearmin = 1995
yearmax = 2019
years = np.arange(1995,2019+1,1)

### Read in obs
latobs,lonobs,levobs,obs = ER.read_ERA5LOWS_monthlyHindcast(variq,'/work/Zachary.Labe/Data/',monthhind,slicenan,'surface')
 
### Calculate time periods
if monthhind == 'Jan':
    if slicemonth == 'FMA':
        obs_time = np.nanmean(obs[:,1:4,:,:],axis=1)
    elif slicemonth == 'MA':
        obs_time = np.nanmean(obs[:,2:4,:,:],axis=1)
    elif slicemonth == 'April':
        obs_time = np.nanmean(obs[:,3:4,:,:],axis=1)
elif monthhind == 'Feb':
    if slicemonth == 'FMA':
        obs_time = np.nanmean(obs[:,0:3,:,:],axis=1)
    elif slicemonth == 'MA':
        obs_time = np.nanmean(obs[:,1:3,:,:],axis=1)
    elif slicemonth == 'April':
        obs_time = np.nanmean(obs[:,2:3,:,:],axis=1)
elif monthhind == 'Mar':
    if slicemonth == 'FMA':
        print(ValueError('FULL TIME PERIOD NOT AVAILABLE!'))
        sys.exit()
    elif slicemonth == 'MA':
        obs_time = np.nanmean(obs[:,0:2,:,:],axis=1)
    elif slicemonth == 'April':
        obs_time = np.nanmean(obs[:,1:2,:,:],axis=1)

### Calculate decadal linear trend
trendobs = UT.linearTrendR(obs_time,years,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################   
### Read in initialized run (low)
lath,lonh,spearlow = SPlow.read_SPEAR_LOW_Hindcast33lev(directorydata,variq,monthhind,slicenan,15)

### Calculate time periods
if monthhind == 'Jan':
    if slicemonth == 'FMA':
        spear_timelow = np.nanmean(spearlow[:,:,1:4,:,:],axis=2)
    elif slicemonth == 'MA':
        spear_timelow = np.nanmean(spearlow[:,:,2:4,:,:],axis=2)
    elif slicemonth == 'April':
        spear_timelow = np.nanmean(spearlow[:,:,3:4,:,:],axis=2)
elif monthhind == 'Feb':
    if slicemonth == 'FMA':
        spear_timelow = np.nanmean(spearlow[:,:,0:3,:,:],axis=2)
    elif slicemonth == 'MA':
        spear_timelow = np.nanmean(spearlow[:,:,1:3,:,:],axis=2)
    elif slicemonth == 'April':
        spear_timelow = np.nanmean(spearlow[:,:,2:3,:,:],axis=2)
elif monthhind == 'Mar':
    if slicemonth == 'FMA':
        print(ValueError('FULL TIME PERIOD NOT AVAILABLE!'))
        sys.exit()
    elif slicemonth == 'MA':
        spear_timelow = np.nanmean(spearlow[:,:,0:2,:,:],axis=2)
    elif slicemonth == 'April':
        spear_timelow = np.nanmean(spearlow[:,:,1:2,:,:],axis=2)
        
### Calculate means across the ensembles
spear_hindlow = np.nanmean(spear_timelow[:,:,:,:],axis=1)
trendhindlow = UT.linearTrendR(spear_hindlow,years,'surface',yearmin,yearmax) * 10.

###############################################################################
###############################################################################
###############################################################################   
### Read in initialized run (high)
lath,lonh,spearhigh = SPhigh.read_SPEAR_LOW_Hindcast65lev(directorydata,variq,monthhind,slicenan,15)

### Calculate time periods
if monthhind == 'Jan':
    if slicemonth == 'FMA':
        spear_timehigh = np.nanmean(spearhigh[:,:,1:4,:,:],axis=2)
    elif slicemonth == 'MA':
        spear_timehigh = np.nanmean(spearhigh[:,:,2:4,:,:],axis=2)
    elif slicemonth == 'April':
        spear_timehigh = np.nanmean(spearhigh[:,:,3:4,:,:],axis=2)
elif monthhind == 'Feb':
    if slicemonth == 'FMA':
        spear_timehigh = np.nanmean(spearhigh[:,:,0:3,:,:],axis=2)
    elif slicemonth == 'MA':
        spear_timehigh = np.nanmean(spearhigh[:,:,1:3,:,:],axis=2)
    elif slicemonth == 'April':
        spear_timehigh = np.nanmean(spearhigh[:,:,2:3,:,:],axis=2)
elif monthhind == 'Mar':
    if slicemonth == 'FMA':
        print(ValueError('FULL TIME PERIOD NOT AVAILABLE!'))
        sys.exit()
    elif slicemonth == 'MA':
        spear_timehigh = np.nanmean(spearhigh[:,:,0:2,:,:],axis=2)
    elif slicemonth == 'April':
        spear_timehigh = np.nanmean(spearhigh[:,:,1:2,:,:],axis=2)
        
### Calculate means across the ensembles
spear_hindhigh = np.nanmean(spear_timehigh[:,:,:,:],axis=1)
trendhindhigh = UT.linearTrendR(spear_hindhigh,years,'surface',yearmin,yearmax) * 10.

###############################################################################
###############################################################################
###############################################################################   
### Read in AMIP 
lat,lon,level,enseall,yearsf = FA.read_FACTS_Experi('spear','spear',variq,slicemonth,4,slicenan,'surface')
spearamipens = UT.linearTrend(enseall,yearsf,'surface',yearmin,yearmax)*10.
spearamip = np.nanmean(spearamipens,axis=0)
###############################################################################
###############################################################################
###############################################################################   
### Read in coupled run
lat,lon,spear = SPM.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',variq,
                                        slicemonth,4,slicenan,30,'all')

trendspear = UT.linearTrend(spear,years,'surface',yearmin,yearmax)*10.
trendspearmean = np.nanmean(trendspear[:,:,:],axis=0)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different AMIPS to compare with SPEAR
fig = plt.figure(figsize=(9,3))

if variq == 'T2M':
    label = r'\textbf{Trend [$^{\circ}$C/decade]}'
    limit = np.arange(-1,1.01,0.02)
    barlim = np.round(np.arange(-1,2,1),2)
elif variq == 'SLP':
    label = r'\textbf{Trend [hPa/decade]}'
    limit = np.arange(-1,1.01,0.02)
    barlim = np.round(np.arange(-1,2,1),2)
elif variq == 'Z500':
    label = r'\textbf{Trend [m/decade]}'
    limit = np.arange(-15,15.01,0.5)
    barlim = np.round(np.arange(-15,16,5),2)
elif variq == 'Z200':
    label = r'\textbf{Trend [m/decade]}'
    limit = np.arange(-15,15.01,0.5)
    barlim = np.round(np.arange(-15,16,5),2)
elif variq == 'Z50':
    label = r'\textbf{Trend [m/decade]}'
    limit = np.arange(-15,15.01,0.5)
    barlim = np.round(np.arange(-15,16,5),2)
elif variq == 'Z30':
    label = r'\textbf{Trend [m/decade]}'
    limit = np.arange(-15,15.01,0.5)
    barlim = np.round(np.arange(-15,16,5),2)
elif variq == 'SST':
    label = r'\textbf{Trend [$^{\circ}$C/decade]}'
    limit = np.arange(-1,1.01,0.02)
    barlim = np.round(np.arange(-1,2,1),2)
elif variq == 'SNOW':
    label = r'\textbf{Trend [m/decade]}'    
    limit = np.arange(-0.03,0.031,0.001)
    barlim = np.round(np.arange(-0.03,0.031,0.03),2)

plotdata = [trendobs,trendhindhigh,trendhindlow,spearamip,trendspearmean]
plotlat = [latobs,lath,lath,lat,lat]
plotlon = [lonobs,lonh,lonh,lon,lon]

for i in range(len(plotdata)):
    ax = plt.subplot(1,5,i+1)
    
    var = plotdata[i]
    lat1 = plotlat[i]
    lon1 = plotlon[i]
    
    m = Basemap(projection='ortho',lon_0=250,lat_0=90,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='dimgrey',linewidth=1)
    m.drawstates(color='dimgrey',linewidth=0.5)
    m.drawcountries(color='dimgrey',linewidth=0.5)
    
    circle = m.drawmapboundary(fill_color='white',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    var, lons_cyclic = addcyclic(var, lon1)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
    x, y = m(lon2d, lat2d)
    
    cs1 = m.contourf(x,y,var,limit,extend='both')
    
    cs1.set_cmap(cmocean.cm.balance)
    
    ## Box 1
    la1 = 43
    la2 = 60
    lo1 = 240
    lo2 = 295
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='gold', linewidth=2, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='gold', linewidth=2, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=2,color='gold',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=2,color='gold',zorder=4)
    
    plt.title(r'%s' % experi1[i],fontsize=8,color='k')
    ax.annotate(r'\textbf{[%s]}' % (letters[i]),xy=(0,0),xytext=(0.03,0.90),
              textcoords='axes fraction',color='k',fontsize=9,
              rotation=40,ha='center',va='center')
    
fig.suptitle(r'\textbf{%s; initialized on %s}' % (slicemonth,monthhind),fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.305,0.1,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
        
plt.savefig(directoryfigure + 'TrendsMap_%s_SPEAR_init-%s_%s_%s-%s.png'% (slicemonth,monthhind,variq,yearmin,yearmax),dpi=300)
