"""
Compare trends in the data for the AMIPS from 1979 - 2020 for only SPEAR
along with its coupled run
 
Author    : Zachary M. Labe
Date      : 21 November 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools
import read_ERA5_monthly1x1 as ER
import read_SPEAR_MED as SP
import read_FACTS_AMIPS as FA

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/AMIPs/TrendsAMIPS/ColdestMember/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
scenario = ['spear_obs_rf']
model_spear = ['spear']
experi1 = ['ERA5','SPEAR_MED-AMIP','SPEAR_MED']
experi2 = ['SPEAR-AMIP','SPEAR_MED','DIFFERENCE']
model = [model_spear]
exxp = 1
modelunravel = list(itertools.chain(*model))
variq = 'T2M'
slicemonth = 'FMA'
slicenan = 'nan'
trendper = 'all'
if trendper == 'all':
    timeslice = '1979-2020'
elif trendper == 'rec':
    timeslice = '2000-2020'
years = np.arange(1979,2020+1,1)

loadens_amip,min_amip = np.load(directorydata + 'EnsembleMembers_SimilarTrends_AMIPSPEAR_T2M_%s_1979-2020.npy' % slicemonth,
                        allow_pickle=True)
loadens_coupled,min_coupled = np.load(directorydata + 'EnsembleMembers_SimilarTrends_CoupledSPEAR_T2M_%s_1979-2020.npy' % slicemonth,
                        allow_pickle=True)

latobs,lonobs,levobs,var = ER.read_ERA5_monthly1x1(variq,'/work/Zachary.Labe/Data/',slicemonth,
                                                    years,3,True,slicenan,'surface')

if variq == 'SST':
    var[np.where(np.isnan(var))] = 0.

### Calculate decadal linear trend
if trendper == 'all':
    yearmin = 1979
    yearmax = 2020
elif trendper == 'rec':
    yearmin = 2000
    yearmax = 2020   
trendobs = UT.linearTrendR(var,years,'surface',yearmin,yearmax)*10

###############################################################################
###############################################################################
###############################################################################   
### Read in AMIP 
lat,lon,levamip,enseall,yearsf = FA.read_FACTS_Experi('spear_obs_rf','SPEAR',variq,slicemonth,4,slicenan,'surface')
spearamipens = UT.linearTrend(enseall,yearsf,'surface',yearmin,yearmax)*10.
spearamip = spearamipens[min_amip]
###############################################################################
###############################################################################
###############################################################################   
### Read in coupled run
lat1s,lon1s,spear = SP.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',variq,
                                        slicemonth,4,
                                        slicenan,30,'all')

trendspear = UT.linearTrend(spear,years,'surface',yearmin,yearmax)*10.
trendspearmean = trendspear[min_coupled]

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different AMIPS to compare with SPEAR
fig = plt.figure(figsize=(9,3.5))

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

plotdata = [trendobs,spearamip,trendspearmean]
plotlat = [latobs,lat,lat]
plotlon = [lonobs,lon,lon]

for i in range(len(plotdata)):
    ax = plt.subplot(1,3,i+1)
    
    var = plotdata[i]
    lat1 = plotlat[i]
    lon1 = plotlon[i]
    
    m = Basemap(projection='ortho',lon_0=250,lat_0=40,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='dimgrey',linewidth=1)
    m.drawstates(color='dimgrey',linewidth=0.5)
    m.drawcountries(color='dimgrey',linewidth=0.5)
    
    circle = m.drawmapboundary(fill_color='white',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    lon2,lat2 = np.meshgrid(lon1,lat1)
    
    cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
    
    cs1.set_cmap(cmocean.cm.balance)
    
    ## Box 2
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
    
fig.suptitle(r'\textbf{%s; %s}' % (slicemonth,timeslice),fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.305,0.08,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
        
plt.savefig(directoryfigure + 'TrendsMap_%s_AMIPsCoupled_%s_%s_ColdestMember.png' % (slicemonth,variq,trendper),dpi=300)
