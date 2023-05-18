"""
Plot monthly trends in the Arctic

Author    : Zachary M. Labe
Date      : 27 June 022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import palettable.cubehelix as ww
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
from scipy.interpolate import griddata as g
import scipy.stats as sts
import cmasher as cmr

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
directorydata = '/work/Zachary.Labe/Data/' 
directorydataAMIP = '/work/Zachary.Labe/Research/Attribution_SpringNA/Data/v1/AMIPs/'

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p",
            "q","r","s","t","u","v","w","x","y","z"]
modelname = ['ECHAM5','CAM5','SPEAR']
variq = 'T2M'
sliceperiod = 'none'
years = np.arange(1979,2021+1,1)
sliceshape = 4
slicenan = 'nan'
addclimo = True

### Assemble data for plotting
plottype = 'spear'
latpolar = 60

# ### Temporary regridding function
# def regrid(lat11,lon11,lat21,lon21,var,years):
#     """
#     Interpolated on selected grid. Reads ERA5 in as 4d with 
#     [year,month,lat,lon]
#     """
    
#     lon1,lat1 = np.meshgrid(lon11,lat11)
#     lon2,lat2 = np.meshgrid(lon21,lat21)
    
#     varn_re = np.reshape(var,(var.shape[0],(lat1.shape[0]*lon1.shape[1])))   
#     varn = np.empty((var.shape[0],lat2.shape[0],lon2.shape[1]))
    
#     print('Completed: Start regridding process:')
#     for i in range(years.shape[0]):
#         z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,:],(lat2,lon2),method='linear')
#         varn[i,:,:] = z
#         print('Completed: Month %s Regridding---' % (years[i]))
#     return varn

# ### Read data
# latobs,lonobs,levobs,var = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
#                                 years,sliceshape,addclimo,slicenan,'surface')

# ### Calculate decadal linear trend
# yearmin = 1979
# yearmax = 2019
# montrend = np.empty((len(monthq),latobs.shape[0],lonobs.shape[0]))
# for i in range(len(montrend)):
#     montrend[i,:,:] = UT.linearTrendR(var[:,i,:,:],years,'surface',yearmin,yearmax)*10

# ### Read in trend data from AMIPS
# amip = np.load(directorydataAMIP + 'AMIP_AMIPSPEAR_MonthlyData_AllMonths_T2M.npz')

# trend_echam = amip['trend_echam']
# lat_echam = amip['lat1_echam']
# lon_echam = amip['lon1_echam']

# trend_cam = amip['trend_cam']
# lat_cam = amip['lat1_cam']
# lon_cam = amip['lon1_cam']

# trend_spear = amip['trend_spear']
# lat_spear = amip['lat1_spear']
# lon_spear = amip['lon1_spear']

# ### Calculate ensemble means
# mean_echam = np.nanmean(trend_echam,axis=1)
# mean_cam = np.nanmean(trend_cam,axis=1)
# mean_spear = np.nanmean(trend_spear,axis=1)

# ### Regrid to match climate model data
# obs_ech = regrid(latobs,lonobs,lat_echam,lon_echam,montrend,np.arange(12))
# obs_cam = regrid(latobs,lonobs,lat_cam,lon_cam,montrend,np.arange(12))
# obs_spe = regrid(latobs,lonobs,lat_spear,lon_spear,montrend,np.arange(12))

# ### Calculate difference
# diff_echam = mean_echam - obs_ech
# diff_cam = mean_cam - obs_cam
# diff_spear = mean_spear - obs_spe

### Calculate regional averaging
def regAvg(lat,lon,data,region,latpolar):
    if region == 'Arctic':
        latq = np.where(lat >= latpolar)[0]
        latnew = lat[latq]
        regionnew = data[:,latq,:]
        lon2new,lat2new = np.meshgrid(lon,latnew)
        ave = UT.calc_weightedAve(regionnew,lat2new) 
    return ave

ave_obs = regAvg(latobs,lonobs,montrend,'Arctic',latpolar)

ave_obs_echam = regAvg(lat_echam,lon_echam,obs_ech,'Arctic',latpolar)
ave_mod_echam = regAvg(lat_echam,lon_echam,mean_echam,'Arctic',latpolar)

ave_obs_cam = regAvg(lat_cam,lon_cam,obs_cam,'Arctic',latpolar)
ave_mod_cam = regAvg(lat_cam,lon_cam,mean_cam,'Arctic',latpolar)

ave_obs_spear = regAvg(lat_spear,lon_spear,obs_spe,'Arctic',latpolar)
ave_mod_spear = regAvg(lat_spear,lon_spear,mean_spear,'Arctic',latpolar)

### Calculate regional warming for each ensemble member
ave_modens_echam = np.empty((trend_echam.shape[1],len(monthq)))
for i in range(trend_echam.shape[1]):
    ave_modens_echam[i,:] = regAvg(lat_echam,lon_echam,trend_echam[:,i,:,:],'Arctic',latpolar)

ave_modens_cam = np.empty((trend_cam.shape[1],len(monthq)))
for i in range(trend_cam.shape[1]):
    ave_modens_cam[i,:] = regAvg(lat_cam,lon_cam,trend_cam[:,i,:,:],'Arctic',latpolar)
    
ave_modens_spear = np.empty((trend_spear.shape[1],len(monthq)))
for i in range(trend_spear.shape[1]):
    ave_modens_spear[i,:] = regAvg(lat_spear,lon_spear,trend_spear[:,i,:,:],'Arctic',latpolar)
    
min_echam = np.nanmin(ave_modens_echam,axis=0)
max_echam = np.nanmax(ave_modens_echam,axis=0)
min_cam = np.nanmin(ave_modens_cam,axis=0)
max_cam = np.nanmax(ave_modens_cam,axis=0)
min_spear = np.nanmin(ave_modens_spear,axis=0)
max_spear = np.nanmax(ave_modens_spear,axis=0)

means = [ave_mod_echam,ave_mod_cam,ave_mod_spear]
minall = [min_echam,min_cam,min_spear]
maxall = [max_echam,max_cam,max_spear]
    
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ### Plot subplot of monthly trends
# limit = np.arange(-2,2.01,0.1)
# barlim = np.round(np.arange(-2,3,1),2)

# fig = plt.figure(figsize=(8,2.5))
# label = r'\textbf{T2M Trend [$^{\circ}$C/decade] for %s--%s}' % (yearmin,yearmax)

# for i in range(len(monthq)*2):
#     ax = plt.subplot(2,12,i+1)
    
#     if i < 12:
#         if plottype == 'echam':
#             var = obs_ech[i]
#             lat1 = lat_echam
#             lon1 = lon_echam
#             modelnameq = modelname[0]
#         elif plottype == 'cam':
#             var = obs_cam[i]
#             lat1 = lat_cam
#             lon1 = lon_cam
#             modelnameq = modelname[1]
#         elif plottype == 'spear':
#             var = obs_spe[i]
#             lat1 = lat_spear
#             lon1 = lon_spear
#             modelnameq = modelname[2]
#         else:
#             print(ValueError('SOMETHING IS WRONG WITH AMIP NAME!'))
#             sys.exit()
#     elif i >=12:
#         if plottype == 'echam':
#             var = mean_echam[i-12]
#             lat1 = lat_echam
#             lon1 = lon_echam
#             modelnameq = modelname[0]
#         elif plottype == 'cam':
#             var = mean_cam[i-12]
#             lat1 = lat_cam
#             lon1 = lon_cam
#             modelnameq = modelname[1]
#         elif plottype == 'spear':
#             var = mean_spear[i-12]
#             lat1 = lat_spear
#             lon1 = lon_spear
#             modelnameq = modelname[2]
#         else:
#             print(ValueError('SOMETHING IS WRONG WITH AMIP NAME!'))
#             sys.exit()
    
#     m = Basemap(projection='npstere',boundinglat=57,lon_0=270,resolution='l',
#                 area_thresh=10000,round=True)
#     m.drawcoastlines(color='dimgrey',linewidth=0.3)
     
#     if plottype == 'spear':
#         circle = m.drawmapboundary(fill_color='white',color='dimgray',
#                           linewidth=0.7)
#         circle.set_clip_on(False)
        
#         lon2,lat2 = np.meshgrid(lon1,lat1)
#         cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
#     else:
#         var, lons_cyclic = addcyclic(var, lon1)
#         var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
#         lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
#         x, y = m(lon2d, lat2d)
        
#         circle = m.drawmapboundary(fill_color='white',color='dimgray',
#                           linewidth=0.7)
#         circle.set_clip_on(False)
        
#         cs1 = m.contourf(x,y,var,limit,extend='both')        
    
#     cs1.set_cmap(cmocean.cm.balance)
    
#     ax.annotate(r'\textbf{[%s]}' % (letters[i]),xy=(0,0),xytext=(0.03,0.90),
#               textcoords='axes fraction',color='k',fontsize=5,
#               rotation=40,ha='center',va='center')
    
#     if i < 12:
#         plt.title(r'\textbf{%s}' % monthq[i],fontsize=12,color='darkgrey')
#     if i == 0:
#         ax.annotate(r'\textbf{ERA5}',xy=(0,0),xytext=(-0.25,0.50),
#                   textcoords='axes fraction',color='w',fontsize=12,
#                   rotation=90,ha='center',va='center')
#     if i == 12:
#         ax.annotate(r'\textbf{%s}' % modelnameq,xy=(0,0),xytext=(-0.25,0.50),
#                   textcoords='axes fraction',color='w',fontsize=12,
#                   rotation=90,ha='center',va='center')
    
# cbar_ax1 = fig.add_axes([0.36,0.16,0.3,0.03])                
# cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
#                     extend='both',extendfrac=0.07,drawedges=False)
# cbar1.set_label(label,fontsize=7,color='darkgrey',labelpad=1.4)  
# cbar1.set_ticks(barlim)
# cbar1.set_ticklabels(list(map(str,barlim)))
# cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
# cbar1.outline.set_edgecolor('darkgrey')
# plt.subplots_adjust(wspace=0,hspace=-0.2)
# # plt.tight_layout()
    
# plt.savefig(directoryfigure + 'PRESENT_TrendsMap_Monthly-Arctic-%s_%s_%s-%s_%s.png' % (plottype,sliceperiod,yearmin,yearmax,variq),dpi=600)

# ###############################################################################
# ###############################################################################
# ###############################################################################
# limit = np.arange(-1,1.01,0.01)
# barlim = np.round(np.arange(-1,2,1),2)

# fig = plt.figure(figsize=(8,1.5))
# label = r'\textbf{T2M Trend Difference [$^{\circ}$C/decade] for %s--%s}' % (yearmin,yearmax)

# for i in range(len(monthq)):
#     ax = plt.subplot(1,12,i+1)
    
#     if i < 12:
#         if plottype == 'echam':
#             var = diff_echam[i]
#             lat1 = lat_echam
#             lon1 = lon_echam
#             modelnameq = modelname[0]
#         elif plottype == 'cam':
#             var = diff_cam[i]
#             lat1 = lat_cam
#             lon1 = lon_cam
#             modelnameq = modelname[1]
#         elif plottype == 'spear':
#             var = diff_spear[i]
#             lat1 = lat_spear
#             lon1 = lon_spear
#             modelnameq = modelname[2]
#         else:
#             print(ValueError('SOMETHING IS WRONG WITH AMIP NAME!'))
#             sys.exit()
    
#     m = Basemap(projection='npstere',boundinglat=57,lon_0=270,resolution='l',
#                 area_thresh=10000,round=True)
#     m.drawcoastlines(color='dimgrey',linewidth=0.3)
        
#     if plottype == 'spear':
#         circle = m.drawmapboundary(fill_color='white',color='dimgray',
#                           linewidth=0.7)
#         circle.set_clip_on(False)
        
#         lon2,lat2 = np.meshgrid(lon1,lat1)
#         cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
#     else:
#         var, lons_cyclic = addcyclic(var, lon1)
#         var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
#         lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
#         x, y = m(lon2d, lat2d)
        
#         circle = m.drawmapboundary(fill_color='white',color='dimgray',
#                           linewidth=0.7)
#         circle.set_clip_on(False)
        
#         cs1 = m.contourf(x,y,var,limit,extend='both')  
    
#     cs1.set_cmap(cmocean.cm.balance)
    
#     ax.annotate(r'\textbf{[%s]}' % (letters[i]),xy=(0,0),xytext=(0.03,0.90),
#               textcoords='axes fraction',color='k',fontsize=5,
#               rotation=40,ha='center',va='center')
    
#     if i < 12:
#         plt.title(r'\textbf{%s}' % monthq[i],fontsize=12,color='darkgrey')
#     if i == 0:
#         ax.annotate(r'\textbf{%s-DIFF}' % modelnameq,xy=(0,0),xytext=(-0.25,0.50),
#                   textcoords='axes fraction',color='w',fontsize=8,
#                   rotation=90,ha='center',va='center')
    
# cbar_ax1 = fig.add_axes([0.36,0.18,0.3,0.03])                
# cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
#                     extend='both',extendfrac=0.07,drawedges=False)
# cbar1.set_label(label,fontsize=7,color='darkgrey',labelpad=1.4)  
# cbar1.set_ticks(barlim)
# cbar1.set_ticklabels(list(map(str,barlim)))
# cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
# cbar1.outline.set_edgecolor('darkgrey')
# plt.subplots_adjust(wspace=0)
# # plt.tight_layout()
    
# plt.savefig(directoryfigure + 'PRESENT_TrendsMap_Monthly-Arctic-%s_%s_%s-%s_%s-DIFFERENCE.png' % (plottype,sliceperiod,yearmin,yearmax,variq),dpi=600)

###############################################################################
###############################################################################
###############################################################################
### Plot mean trend over polar cap
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([])

fig = plt.figure(figsize=(8,2))
ax = plt.subplot(111)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(0)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey',
               labelbottom='off',bottom='off')
ax.tick_params(axis = "x", which = "both", bottom = False, top = False)

ccc=cmr.apple(np.linspace(0.4,1,len(modelname)))
for i in range(len(monthq)):
    plt.scatter(i,ave_obs[i],s=50,c='aqua',edgecolor='aqua',zorder=5,clip_on=False)
    plt.errorbar(i-0.2,means[0][i],
                  yerr=np.array([[means[0][i]-minall[0][i],maxall[0][i]-means[0][i]]]).T,
                  color=ccc[0],linewidth=0.5,capthick=2,capsize=4,clip_on=False)
    plt.errorbar(i,means[1][i],
                  yerr=np.array([[means[1][i]-minall[1][i],maxall[1][i]-means[1][i]]]).T,
                  color=ccc[1],linewidth=0.5,capthick=2,capsize=4,clip_on=False)
    plt.errorbar(i+0.2,means[2][i],
                  yerr=np.array([[means[2][i]-minall[2][i],maxall[2][i]-means[2][i]]]).T,
                  color=ccc[2],linewidth=0.5,capthick=2,capsize=4,clip_on=False)
    
plt.text(3,1.5,r'\textbf{ECHAM5}',fontsize=10,color=ccc[0])
plt.text(4.3,1.5,r'\textbf{CAM5}',fontsize=10,color=ccc[1])
plt.text(5.3,1.5,r'\textbf{SPEAR}',fontsize=10,color=ccc[2])
plt.text(6.3,1.5,r'\textbf{ERA5 [obs]}',fontsize=10,color='aqua')

plt.ylabel(r'\textbf{T2M Trend - Arctic - [$^{\circ}$C/decade]}',color='w',fontsize=5)    
plt.xticks(np.arange(0,12,1),monthq,fontsize=7,color='w')
plt.yticks(np.arange(-5,5.1,0.25),map(str,np.round(np.arange(-5,5.1,0.25),2)),fontsize=7)
plt.xlim([-1,12])
plt.ylim([0,1.5])
plt.tight_layout()
plt.savefig(directoryfigure + 'PRESENT_meanTrends_AMIPS_Arctic-%s_%s_%s_%s.png' % (latpolar,variq,yearmin,yearmax),dpi=600)


