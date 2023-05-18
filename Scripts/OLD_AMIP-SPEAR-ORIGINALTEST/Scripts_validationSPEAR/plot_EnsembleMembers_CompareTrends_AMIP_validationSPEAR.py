"""
Compare ensemble members with similar trends or even colder
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

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/validationSPEAR/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/'

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','spear']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
model = [model_1880s_rf,model_obs_rf,model_spear]
model_spear = ['spear']
slicemonth = 'FMA'
variq = 'T2M'
trendper= 'all'
modelunravel = list(itertools.chain(*model))
yearsobs_all = np.arange(1979,2019+1,1)
yearsobs_rec = np.arange(2000,2019+1,1)

### Read in data
datatrends = np.load(directorydata + 'MapsOfTrends_SPEARval_AMIPs_%s_%s_%s.npz' % (slicemonth,
                                                            variq,trendper),allow_pickle=True)
ens_gridold = datatrends['trends']
lat = datatrends['lat']
lon = datatrends['lon']

### Prepare individual amips to compare
echam = np.asarray(ens_gridold[1][0])
cam = np.asarray(ens_gridold[1][1])
spear = np.asarray(ens_gridold[2][0])
alllats = lat[1]
lat_ech = np.asarray(alllats[0])
lat_cam = np.asarray(alllats[1])
lat_spe = np.asarray(lat[2][0])
alllons = lon[1]
lon_ech = np.asarray(alllons[0])
lon_cam = np.asarray(alllons[1])
lon_spe = np.asarray(lon[2][0])

### Read in ensemble members
ensq = np.load(directorydata + 'EnsembleMembers_SimilarTrends_validationSPEAR_T2M_%s_%s.npz' % (slicemonth,trendper),allow_pickle=True)
echam_q = ensq['echam']
cam_q = ensq['cam']
spear_q = ensq['spear']

### Only ensembles with trends similar
echam_only = echam[echam_q,:,:]
cam_only = cam[cam_q,:,:]
spear_only = spear[spear_q,:,:]

###############################################################################
###############################################################################
###############################################################################
### Plot subplot ECHAM5
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

fig = plt.figure(figsize=(10,3))

for i in range(len(echam_only)):
    ax = plt.subplot(1,len(echam_only),i+1)
    
    var = echam_only[i]
    lat1 = lat_ech
    lon1 = lon_ech
    
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
    la2 = 60
    lo1 = 240
    lo2 = 295
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1,color='gold',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1,color='gold',zorder=4)
    
    plt.title(r'Ensemble Member - %s' % (echam_q[i]+1),fontsize=5,color='k')
    
fig.suptitle(r'\textbf{%s - ECHAM5(%s)}' % (slicemonth,len(echam)),fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.33,0.1,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
        
plt.savefig(directoryfigure + 'EnsembleMembers_CompareTrends_ECHAM5_%s_AMIPs_%s_%s.png'% (slicemonth,variq,trendper),dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot CAM5
limit = np.arange(-1,1.01,0.02)
barlim = np.round(np.arange(-1,2,1),2)
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

fig = plt.figure(figsize=(10,5))

for i in range(len(cam_only)):
    ax = plt.subplot(1,len(cam_only),i+1)
    
    var = cam_only[i]
    lat1 = lat_cam
    lon1 = lon_cam
    
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
    la2 = 60
    lo1 = 240
    lo2 = 295
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1,color='gold',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1,color='gold',zorder=4)
    
    plt.title(r'Ensemble Member - %s' % (cam_q[i]+1),fontsize=8,color='k')
    
fig.suptitle(r'\textbf{%s - ESRL-CAM5(%s)}' % (slicemonth,len(cam)),fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.31,0.09,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
        
plt.savefig(directoryfigure + 'EnsembleMembers_CompareTrends_ESRL-CAM5_%s_AMIPs_%s_%s.png'% (slicemonth,variq,trendper),dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot subplot SPEAR
limit = np.arange(-1,1.01,0.02)
barlim = np.round(np.arange(-1,2,1),2)
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

fig = plt.figure(figsize=(8,6))

for i in range(len(spear_only)):
    ax = plt.subplot(1,len(spear_only),i+1)
    
    var = spear_only[i]
    lat1 = lat_spe
    lon1 = lon_spe
    
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
    la2 = 60
    lo1 = 240
    lo2 = 295
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='gold', linewidth=1, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1,color='gold',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1,color='gold',zorder=4)
    
    plt.title(r'Ensemble Member - %s' % (spear_q[i]+1),fontsize=8,color='k')
    
fig.suptitle(r'\textbf{%s - SPEAR(%s)}' % (slicemonth,len(spear)),fontsize=18,color='k')
    
cbar_ax1 = fig.add_axes([0.33,0.1,0.4,0.03])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
        
plt.savefig(directoryfigure + 'EnsembleMembers_CompareTrends_SPEAR_%s_AMIPs_%s_%s.png' % (slicemonth,variq,trendper),dpi=300)
