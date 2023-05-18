"""
Plot trends in FMA from 1979 to 2014 for the AMIP simulations
Author    : Zachary M. Labe
Date      : 6 June 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','amip_clim_polar','spear']
experi = ['amip_1880s_rf','amip_1880s_rf','amip_1880s_rf',
          'amip_obs_rf','amip_obs_rf','amip_obs_rf',
          'amip_clim_polar','amip_clim_polar','amip_clim_polar','amip_clim_polar',
          'spear']
model_1880s_rf = ['CAM4','ECHAM5','ESRL-CAM5']
model_obs_rf = ['CAM4','ECHAM5','ESRL-CAM5']
model_clim_polar = ['CAM4','ECHAM5','ESRL-CAM5','ESRL-GFS']
model_spear = ['spear']
model = [model_1880s_rf,model_obs_rf,model_clim_polar,model_spear]
modelunravel = list(itertools.chain(*model))
variq = 'T2M'
slicemonth= 'April'
trendper = 'recent'
slicenan = 'nan'
datareader = False
calctrends = True

### Read data
if calctrends == True:
    fma = np.load(directorydata + 'AMIP_AllData_%s_%s.npz' % (slicemonth,variq),
                  allow_pickle = True)
    year = np.load(directorydata + 'AMIP_Years_%s_%s.npy' % (slicemonth,variq),
                        allow_pickle=True)
    ens = np.load(directorydata + 'AMIP_EnsembleMembers_%s_%s.npy' % (slicemonth,variq),
                        allow_pickle=True)
    lenyear = np.load(directorydata + 'AMIP_LengthYears_%s_%s.npy' % (slicemonth,variq),
                        allow_pickle=True)
    lat = fma['lat']
    lon = fma['lon']
    amip = fma['amip']
    
    ### Calculate map of linear trend for each ensemble member
    treeall = []
    for sc in range(len(scenario)):
        tree = []
        for ee in range(len(model[sc])):
            ### Calculate decadal linear trend
            yearq = year[sc][ee]
            if trendper == 'all':
                yearmin = 1979
                yearmax = 2014
            elif trendper == 'recent':
                yearmin = 2000
                yearmax = 2014
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
        
    np.savez(directorydata + 'MapsOfTrends_AMIPs_%s_%s_%s.npz' % (slicemonth,
                                                               variq,
                                                               trendper),
             lat=lat,lon=lon,trends=treeall)
    np.savez(directorydata + 'MapsOfTrends_AMIPs_EnsembleMean_%s_%s_%s.npz' % (slicemonth,
                                                               variq,
                                                               trendper),
             lat=lat,lon=lon,trends=enseall)
    
###############################################################################
###############################################################################
###############################################################################    
if datareader == True:
    datatrends = np.load(directorydata + 'MapsOfTrends_AMIPs_EnsembleMean_%s_%s_%s.npz' % (slicemonth,
                                                               variq,
                                                               trendper),
                         allow_pickle=True)
    enseall = datatrends['trends']
    lat = datatrends['lat']
    lon = datatrends['lon']

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-1,1.01,0.02)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(10,10))
label = r'\textbf{Trend [$^{\circ}$C/decade]}'

plotdata = [enseall[0][0],enseall[0][1],enseall[0][2],
            enseall[1][0],enseall[1][1],enseall[1][2],
            enseall[2][0],enseall[2][1],enseall[2][2],enseall[2][3],
            enseall[3][0]]
plotlat = [lat[0][0],lat[0][1],lat[0][2],
            lat[1][0],lat[1][1],lat[1][2],
            lat[2][0],lat[2][1],lat[2][2],lat[2][3],
            lat[3][0]]
plotlon = [lon[0][0],lon[0][1],lon[0][2],
            lon[1][0],lon[1][1],lon[1][2],
            lon[2][0],lon[2][1],lon[2][2],lon[2][3],
            lon[3][0]]
plotindex = [1,2,3,5,6,7,9,10,11,12,13]

for i in range(len(plotdata)):
    ax = plt.subplot(4,4,plotindex[i])
    
    var = plotdata[i]
    lat1 = plotlat[i]
    lon1 = plotlon[i]
    
    m = Basemap(projection='ortho',lon_0=250,lat_0=40,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='dimgrey',linewidth=1)
    m.drawstates(color='dimgrey',linewidth=0.5)
    m.drawcountries(color='dimgrey',linewidth=0.5)
        
    # var, lons_cyclic = addcyclic(var, lon1)
    # var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    # lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
    # x, y = m(lon2d, lat2d)
    
    circle = m.drawmapboundary(fill_color='white',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    lon2,lat2 = np.meshgrid(lon1,lat1)
    
    cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
    
    cs1.set_cmap(cmocean.cm.balance)
    
    ## Box 2
    la1 = 43
    la2 = 63
    lo1 = 238
    lo2 = 270
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1,color='gold',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1,color='gold',zorder=4)
    
    plt.title(r'%s' % modelunravel[i],fontsize=8,color='k')
    
    if any([i == 0,i == 3,i == 6,i == 10]):
        ax.annotate(r'\textbf{%s}' % experi[i],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=17,rotation=90,ha='center',va='center')
    
fig.suptitle(r'\textbf{APRIL}',fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.4,0.05,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
        
plt.savefig(directoryfigure + 'TrendsMap_%s_AMIPs_%s_%sz.png'% (slicemonth,variq,trendper),dpi=300)
