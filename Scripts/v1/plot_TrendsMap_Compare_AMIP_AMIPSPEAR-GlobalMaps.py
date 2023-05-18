"""
Compare trends in the data for the AMIPS from 1979 - 2019 (global maps)
 
Author    : Zachary M. Labe
Date      : 20 June 2022
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

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/AMIPS/AMIPtrends/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','spear_obs_rf']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
experi = ['ECHAM5','ESRL-CAM5','SPEAR','ERA5']
model = [model_1880s_rf,model_obs_rf,model_spear]
exxp = 1
modelunravel = list(itertools.chain(*model))
variq = 'T2M'
slicemonth = 'annual'
trendper = 'all'

if trendper == 'all':
    timeslice = '1979-2019'
elif trendper == 'rec':
    timeslice = '2000-2019'
slicenan = 'nan'
datareader = True
calctrends = False
years = np.arange(1979,2019+1,1)

latobs,lonobs,levobs,var = ER.read_ERA5_monthly1x1(variq,'/work/Zachary.Labe/Data/',slicemonth,
                                years,3,True,slicenan,'surface')

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2019
trendobs = UT.linearTrendR(var,years,'surface',yearmin,yearmax)*10

### Read data
if calctrends == True:
    fma = np.load(directorydata + 'AMIP_AMIPSPEAR_AllData_%s_%s.npz' % (slicemonth,variq),
                  allow_pickle = True)
    lat = fma['lat']
    lon = fma['lon']
    amip = fma['amip']
    
    ### Calculate map of linear trend for each ensemble member
    treeall = []
    for sc in range(len(scenario)):
        tree = []
        for ee in range(len(model[sc])):
            ### Calculate decadal linear trend
            yearq = years
            if trendper == 'all':
                yearmin = 1979
                yearmax = 2019
            elif trendper == 'recent':
                yearmin = 2000
                yearmax = 2019
            trend = UT.linearTrend(amip[sc][ee],yearq,'surface',
                                   yearmin,yearmax)*10
            tree.append(trend)
            print('Finished for model ---> %s' % model[sc][ee])
        treeall.append(tree)
        print('Finished for experiment ---> %s' % scenario[sc])
            
    ### Calculate ensemble means
    enseall = []
    for sc in range(len(scenario)):
        ense = []
        for ee in range(len(model[sc])):
            ensmean = np.nanmean(treeall[sc][ee],axis=0)
            ense.append(ensmean)
        enseall.append(ense)
        
    np.savez(directorydata + 'MapsOfTrends_AMIPSPEAR_AMIPs_%s_%s_%s.npz' % (slicemonth,
                                                               variq,
                                                               trendper),
             lat=lat,lon=lon,trends=treeall)
    np.savez(directorydata + 'MapsOfTrends_AMIPSPEAR_AMIPs_EnsembleMean_%s_%s_%s.npz' % (slicemonth,
                                                               variq,
                                                               trendper),
             lat=lat,lon=lon,trends=enseall)

###############################################################################
###############################################################################
###############################################################################    
if datareader == True:
    datatrends = np.load(directorydata + 'MapsOfTrends_AMIPSPEAR_AMIPs_EnsembleMean_%s_%s_%s.npz' % (slicemonth,
                                                               variq,
                                                               trendper),
                         allow_pickle=True)
    enseall = datatrends['trends']
    lat = datatrends['lat']
    lon = datatrends['lon']

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
    label = r'\textbf{Trend [$^{\circ}$hPa/decade]}'
    limit = np.arange(-1,1.01,0.02)
    barlim = np.round(np.arange(-1,2,1),2)

plotdata = [enseall[exxp][0],enseall[exxp][1],enseall[-1][0],trendobs]
plotlat = [lat[exxp][0],lat[exxp][1],lat[-1][0],latobs]
plotlon = [lon[exxp][0],lon[exxp][1],lon[-1][0],lonobs]


for i in range(len(experi)):
    ax = plt.subplot(1,4,i+1)
    
    var = plotdata[i]
    lat1 = plotlat[i]
    lon1 = plotlon[i]
    
    m = Basemap(projection='robin',lon_0=0,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='darkgrey',linewidth=0.5)
    
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
    m.plot(lonsslice, latsslice, color='aqua', linewidth=1, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='aqua', linewidth=1, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1,color='aqua',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1,color='aqua',zorder=4)
    
    plt.title(r'%s' % experi[i],fontsize=8,color='k')
    
    if any([i == 0]):
        ax.annotate(r'\textbf{%s}' % scenario[exxp],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=17,rotation=90,ha='center',va='center')
    
fig.suptitle(r'\textbf{%s; %s}' % (slicemonth,timeslice),fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.318,0.1,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
        
plt.savefig(directoryfigure + 'TrendsMap_%s_AMIPs_%s_%s_%s_%s-%s_GlobalMaps.png'% (slicemonth,variq,trendper,scenario[exxp],yearmin,yearmax),dpi=300)
