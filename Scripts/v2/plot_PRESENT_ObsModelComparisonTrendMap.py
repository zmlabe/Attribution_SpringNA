"""
Compare trends in the data for the AMIPS from 1979 - 2020 for only SPEAR
along with its coupled run and FLOR for AGU
 
Author    : Zachary M. Labe
Date      : 8 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools
import read_ERA5_monthlyMEDS as ER
import read_SPEAR_MED as SP
import read_FLOR as FL
import read_FACTS_AMIPS as FA

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

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
scenario = ['spear_obs_rf']
model_spear = ['spear']
experi1 = ['ERA5','SPEAR_MED_AMIP','SPEAR_MED_LE','FLOR_LE']
variq = 'SLP'
slicemonth = 'FMA'
slicenan = 'nan'
trendper = 'all'
if trendper == 'all':
    timeslice = '1979-2020'
elif trendper == 'rec':
    timeslice = '2000-2020'
years = np.arange(1979,2020+1,1)

latobs,lonobs,levobs,var = ER.read_ERA5_monthlyMEDS(variq,'/work/Zachary.Labe/Data/',slicemonth,
                                                    years,3,True,slicenan,'surface')

### Calculate decadal linear trend
if trendper == 'all':
    yearmin = 1979
    yearmax = 2020
trendobs = UT.linearTrendR(var,years,'surface',yearmin,yearmax)*10
###############################################################################
###############################################################################
###############################################################################   
### Read in AMIP 
lat,lon,levamip,enseall,yearsf = FA.read_FACTS_Experi('spear_obs_rf','SPEAR',variq,slicemonth,4,slicenan,'surface')
spearamipens = UT.linearTrend(enseall,yearsf,'surface',yearmin,yearmax)*10.
spearamip = np.nanmean(spearamipens,axis=0)

###############################################################################
###############################################################################
###############################################################################   
### Read in coupled run
lat1s,lon1s,spear = SP.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',variq,
                                        slicemonth,4,
                                        slicenan,30,'all')

trendspear = UT.linearTrend(spear,years,'surface',yearmin,yearmax)*10.
trendspearmean = np.nanmean(trendspear[:,:,:],axis=0)

###############################################################################
###############################################################################
###############################################################################   
### Read in FLOR
lat1s,lon1s,flor = FL.read_FLOR('/work/Zachary.Labe/Data/',variq,
                                        slicemonth,4,
                                        slicenan,30,'satellite')

trendflor = UT.linearTrend(flor,years,'surface',yearmin,yearmax)*10.
trendflormean = np.nanmean(trendflor[:,:,:],axis=0)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different AMIPS to compare with SPEAR
fig = plt.figure(figsize=(9,2.7))

if variq == 'T2M':
    label = r'\textbf{T2M Trends [$^{\circ}$C/Decade] for April %s-%s' % (1979,2020)
    limit = np.arange(-1,1.01,0.01)
    barlim = np.round(np.arange(-1,2,1),2)
elif variq == 'SLP':
    label = r'\textbf{SLP Trends [hPa/Decade] for April %s-%s' % (1979,2020)
    limit = np.arange(-1,1.01,0.01)
    barlim = np.round(np.arange(-1,2,1),2)
elif variq == 'Z500':
    label = r'\textbf{Z500 Trends [m/Decade] for April %s-%s' % (1979,2020)
    limit = np.arange(-20,20.01,0.1)
    barlim = np.round(np.arange(-20,21,10),2)

plotdata = [trendobs,spearamip,trendspearmean,trendflormean]
plotlat = [lat,lat,lat,lat]
plotlon = [lon,lon,lon,lon]

for i in range(len(plotdata)):
    ax = plt.subplot(1,4,i+1)
    
    var = plotdata[i]
    lat1 = plotlat[i]
    lon1 = plotlon[i]
    
    m = Basemap(projection='robin',lon_0=0,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='dimgrey',linewidth=0.5,zorder=12)
    
    circle = m.drawmapboundary(fill_color='k',color='k',linewidth=1)
    circle.set_clip_on(False)
    
    lon2,lat2 = np.meshgrid(lon1,lat1)
    
    cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
    
    cs1.set_cmap(cmocean.cm.balance)
    
    plt.title(r'\textbf{%s}' % experi1[i],fontsize=15,color='w')
    
    # m.drawlsmask(land_color=(0,0,0,0),ocean_color='k',lakes=False,zorder=11)
    
cbar_ax1 = fig.add_axes([0.305,0.2,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='w',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=10)
cbar1.outline.set_edgecolor('darkgrey')
plt.tight_layout()
        
plt.savefig(directoryfigure + 'POSTER_TrendsMap_%s_AMIPsCoupled_%s_%s.png'% (slicemonth,variq,trendper),dpi=600)
