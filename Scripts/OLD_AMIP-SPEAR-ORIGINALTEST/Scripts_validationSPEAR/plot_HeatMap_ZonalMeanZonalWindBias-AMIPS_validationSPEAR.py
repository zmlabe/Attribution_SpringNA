"""
Plot monthly bias of zonal-mean zonal wind in the AMIP experiments

Author    : Zachary M. Labe
Date      : 8 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts
import cmocean
import read_ERA5_monthly1x1 as ER

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/validationSPEAR/' 
directorydata = '/work/Zachary.Labe/Data/' 
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/'

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','spear']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
model = [model_1880s_rf,model_obs_rf,model_spear]
variq = 'U'
level = 'vertical'
slicemonth= 'none'
slicenan = 'nan'
addclimo = True
yearAllData = np.arange(1979,2019+1,1)
months = np.arange(len(monthq))

### Calculate polar cap height
def calc_regionalAve(regionn,lat1,lon1,data):
    if regionn == 'CANMID':
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
    elif regionn == 'PCH':
        la1 = 65
        la2 = 90
        lo1 = 0
        lo2 = 360
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
        
    if data.ndim == 3:
        meanlat = data[:,lat1q,:]
        meanbox = meanlat[:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    elif data.ndim == 4:
        meanlat = data[:,:,lat1q,:]
        meanbox = meanlat[:,:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    elif data.ndim == 5:
        meanlat = data[:,:,:,lat1q,:]
        meanbox = meanlat[:,:,:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    elif data.ndim == 6:
        meanlat = data[:,:,:,:,lat1q,:]
        meanbox = meanlat[:,:,:,:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
        
    ### Calculate timeseries
    mean = UT.calc_weightedAve(meanbox,lat2q)

    return mean

### Read in observational data
latobs,lonobs,levobs,varobs = ER.read_ERA5_monthly1x1(variq,directorydata,slicemonth,
                                                      yearAllData,5,addclimo,slicenan,level)

obspch = calc_regionalAve('PCH',latobs,lonobs,varobs[:-2,:,:,:,:])
obspchf = np.flip(obspch,axis=2)
levobsf = np.flip(levobs)

###############################################################################
###############################################################################
###############################################################################
### Read in bias of zonal-mean zonal wind
biasfile = np.load(directoryoutput + 'AMIP_SPEARval_MonthlyData_U-PCH_Bias.npz',allow_pickle=True)
bias = biasfile['bias']
levels = biasfile['levels']

echam_1880 = bias[0][0]
echam_1880_mean = np.nanmean(echam_1880,axis=(0,1)).transpose()
cam_1880 = bias[0][1]
cam_1880_mean = np.nanmean(cam_1880,axis=(0,1)).transpose()

echam_obs = bias[1][0]
echam_obs_mean = np.nanmean(echam_obs,axis=(0,1)).transpose()
cam_obs = bias[1][1]
cam_obs_mean = np.nanmean(cam_obs,axis=(0,1)).transpose()
spear_obs = bias[2][0]
spear_obs_mean = np.nanmean(spear_obs,axis=(0,1)).transpose()

###############################################################################
###############################################################################
###############################################################################
### Read in zonal-mean zonal wind
datafile = np.load(directoryoutput + 'AMIP_SPEARval_MonthlyData_U-PCH.npz',allow_pickle=True)
u = datafile['u']

echam_1880u = u[0][0]
echam_1880u_mean = np.nanmean(echam_1880u,axis=(0,1)).transpose()
cam_1880u = u[0][1]
cam_1880u_mean = np.nanmean(cam_1880u,axis=(0,1)).transpose()

echam_obsu = u[1][0]
echam_obsu_mean = np.nanmean(echam_obsu,axis=(0,1)).transpose()
cam_obsu = u[1][1]
cam_obsu_mean = np.nanmean(cam_obsu,axis=(0,1)).transpose()
spear_obsu = u[2][0]
spear_obsu_mean = np.nanmean(spear_obsu,axis=(0,1)).transpose()

###############################################################################
###############################################################################
###############################################################################
### Ready plotting data
month2,lev2 = np.meshgrid(months,levels)

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(131)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom=True)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft='on')

cs = plt.pcolormesh(echam_obs_mean,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-10,vmax=10)

cs.set_cmap(cmocean.cm.balance)

plt.yticks(np.arange(0.5,len(levels)+0.5,1),levels.astype(int),ha='right',color='dimgrey',
            va='center',fontsize=6)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)
plt.ylim([2,17])
plt.ylabel(r'\textbf{Pressure (hPa)}',color='k',fontsize=7)

plt.xticks(np.arange(0.5,len(monthq)+0.5,1),monthq,ha='center',color='dimgrey',
            va='center',rotation=90,fontsize=6)
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,12])

plt.title(r'\textbf{ECHAM5-AMIP-OBS}',color='dimgrey',fontsize=11)

###############################################################################
###############################################################################
###############################################################################
ax = plt.subplot(132)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom=True)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft='on')

cs = plt.pcolormesh(cam_obs_mean,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-10,vmax=10)

cs.set_cmap(cmocean.cm.balance)

plt.yticks(np.arange(0.5,len(levels)+0.5,1),levels.astype(int),ha='right',color='dimgrey',
            va='center',fontsize=6)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)
plt.ylim([2,17])

plt.xticks(np.arange(0.5,len(monthq)+0.5,1),monthq,ha='center',color='dimgrey',
            va='center',rotation=90,fontsize=6)
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,12])

plt.title(r'\textbf{CAM5-AMIP-OBS}',color='dimgrey',fontsize=11)

###############################################################################
###############################################################################
###############################################################################
ax = plt.subplot(133)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom=True)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft='on')

cs = plt.pcolormesh(spear_obs_mean,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-10,vmax=10)

cs.set_cmap(cmocean.cm.balance)

plt.yticks(np.arange(0.5,len(levels)+0.5,1),levels.astype(int),ha='right',color='dimgrey',
            va='center',fontsize=6)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)
plt.ylim([2,17])

plt.xticks(np.arange(0.5,len(monthq)+0.5,1),monthq,ha='center',color='dimgrey',
            va='center',rotation=90,fontsize=6)
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,12])

plt.title(r'\textbf{SPEAR-AMIP-OBS}',color='dimgrey',fontsize=11)

fig.suptitle(r'\textbf{BIAS IN ZONAL-MEAN, ZONAL WIND}',color='k')

cbar_ax1 = fig.add_axes([0.318,0.08,0.4,0.03])                
cbar1 = fig.colorbar(cs,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(r'\textbf{[U] m/s}',fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(np.arange(-10,11,2))
cbar1.set_ticklabels(list(map(str,np.arange(-10,11,2))))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')

plt.tight_layout()
plt.subplots_adjust(bottom=0.2)
plt.savefig(directoryfigure + 'HeatMap_ZonalMeanZonalWindBias-AMIPS-obs_validationSPEAR.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(121)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom=True)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft='on')

cs = plt.pcolormesh(echam_1880_mean,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-10,vmax=10)

cs.set_cmap(cmocean.cm.balance)

plt.yticks(np.arange(0.5,len(levels)+0.5,1),levels.astype(int),ha='right',color='dimgrey',
            va='center',fontsize=6)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)
plt.ylim([2,17])
plt.ylabel(r'\textbf{Pressure (hPa)}',color='k',fontsize=7)

plt.xticks(np.arange(0.5,len(monthq)+0.5,1),monthq,ha='center',color='dimgrey',
            va='center',rotation=90,fontsize=6)
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,12])

plt.title(r'\textbf{ECHAM5-AMIP-1880}',color='dimgrey',fontsize=11)

###############################################################################
###############################################################################
###############################################################################
ax = plt.subplot(122)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom=True)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft='on')

cs = plt.pcolormesh(cam_1880_mean,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-10,vmax=10)

cs.set_cmap(cmocean.cm.balance)

plt.ylabel(r'\textbf{Pressure (hPa)}',color='k',fontsize=7)
plt.yticks(np.arange(0.5,len(levels)+0.5,1),levels.astype(int),ha='right',color='dimgrey',
            va='center',fontsize=6)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)
plt.ylim([2,17])

plt.xticks(np.arange(0.5,len(monthq)+0.5,1),monthq,ha='center',color='dimgrey',
            va='center',rotation=90,fontsize=6)
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,12])

plt.title(r'\textbf{CAM5-AMIP-1880}',color='dimgrey',fontsize=11)

cbar_ax1 = fig.add_axes([0.318,0.08,0.4,0.03])                
cbar1 = fig.colorbar(cs,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(r'\textbf{[U] m/s}',fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(np.arange(-10,11,2))
cbar1.set_ticklabels(list(map(str,np.arange(-10,11,2))))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')

plt.tight_layout()
plt.subplots_adjust(bottom=0.2)
plt.savefig(directoryfigure + 'HeatMap_ZonalMeanZonalWindBias-AMIPS-1880_validationSPEAR.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Strongest zonal winds
sortedobs = np.sort(obspchf,axis=0)
strength = np.nanmean(sortedobs[-5:,:,:] - sortedobs[:5,:,:],axis=0)

### Loop through echam
strength_echam = np.empty((len(echam_obsu),months.shape[0],levels.shape[0]))
for i in range(len(strength_echam)):
    sorted_echam = np.sort(echam_obsu[i,:,:,:],axis=0)
    strength_echam[i] = np.nanmean(sorted_echam[-5:,:,:] - sorted_echam[:5,:,:],axis=0)
strength_echamm = np.nanmean(strength_echam,axis=0)
diff_strength_echam = strength_echamm - strength

### Loop through cam
strength_cam = np.empty((len(cam_obsu),months.shape[0],levels.shape[0]))
for i in range(len(strength_cam)):
    sorted_cam = np.sort(cam_obsu[i,:,:,:],axis=0)
    strength_cam[i] = np.nanmean(sorted_cam[-5:,:,:] - sorted_cam[:5,:,:],axis=0)
strength_camm = np.nanmean(strength_cam,axis=0)
diff_strength_cam = strength_camm - strength

### Loop through spear
strength_spear = np.empty((len(spear_obsu),months.shape[0],levels.shape[0]))
for i in range(len(strength_spear)):
    sorted_spear = np.sort(spear_obsu[i,:,:,:],axis=0)
    strength_spear[i] = np.nanmean(sorted_spear[-5:,:,:] - sorted_spear[:5,:,:],axis=0)
strength_spearm = np.nanmean(strength_spear,axis=0)
diff_strength_spear = strength_spearm - strength
                   
###############################################################################
###############################################################################
###############################################################################      
assembleForPlotting = [strength_echamm,strength,diff_strength_echam,
                        strength_camm,strength,diff_strength_cam,
                        strength_spearm,strength,diff_strength_spear]
titlelabels = ['AMIP-OBS','ERA5','DIFFERENCE']
labelmodels = np.repeat(['ECHAM5','CAM5','SPEAR'],3)

###########################################################################
###########################################################################
###########################################################################
##### Plot profiles
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
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
        
### Set limits for contours and colorbars6
limit = np.arange(-30,30.1,0.5)
barlim = np.arange(-30,31,5)
cmap = cmocean.cm.balance
label = r'\textbf{[U] m/s}'
zscale = np.array([1000,700,500,300,200,100,50,30,10])
        
fig = plt.figure(figsize=(10,4))

### Create plot
for p in range(len(assembleForPlotting)):
    ax1 = plt.subplot(3,3,p+1)
    ax1.spines['top'].set_color('dimgrey')
    ax1.spines['right'].set_color('dimgrey')
    ax1.spines['bottom'].set_color('dimgrey')
    ax1.spines['left'].set_color('dimgrey')
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['right'].set_linewidth(2)
    ax1.spines['top'].set_linewidth(2)
    ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='dimgrey')   
    ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='dimgrey')
    plt.gca().axes.get_xaxis().set_visible(True)
    plt.gca().axes.get_yaxis().set_visible(True)
    
    if any([p==0,p==3,p==6]):
        plt.ylabel(r'\textbf{Pressure [hPa]}',color='k',fontsize=7)
    
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    
    var = assembleForPlotting[p]
    
    ### Plot contours
    cs = plt.contourf(month2,lev2,var.transpose(),limit,extend='both')
    # cs1 = plt.contour(timeq,levq,climowinter.transpose(),np.arange(-10,35,2),
    #                   linewidths=2,colors='dimgrey')
    
    cs.set_cmap(cmap)
    
    plt.gca().invert_yaxis()
    plt.yscale('log')
    
    if any([p==6,p==7,p==8]):
        plt.xticks(np.arange(0,11+1,1),monthq,fontsize=4)
    else:
        plt.xticks([])
        
    if p < 3:
        plt.title(r'\textbf{%s}' % titlelabels[p],color='k',fontsize=12)
    
    if any([p==0,p==3,p==6]):
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    else:
        plt.yticks([])
        
    if any([p==2,p==5,p==8]):
        plt.text(11.6,100,r'\textbf{%s}' % labelmodels[p],fontsize=15,rotation=270,
                  color='dimgrey',ha='center',va='center')
    
    plt.xlim([0,11])
    plt.ylim([1000,10])
    plt.minorticks_off()

###########################################################################
plt.tight_layout()
cbar_ax = fig.add_axes([0.32,0.08,0.4,0.03])                    
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(label,fontsize=11,color='dimgrey',labelpad=1.4)  

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001,labelsize=7)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(bottom=0.17,hspace=0.12,wspace=0.08)    
plt.savefig(directoryfigure + 'HeatMap_STRONGWEAK-PV_ZonalMeanZonalWindBias-AMIPS-obs_validationSPEAR.png',dpi=300)
print('Completed: Script done!')

###############################################################################
###############################################################################
###############################################################################
### Strongest/weakest PV based on U50
lev50 = np.where((levels == 50))[0]
obs50 = obspchf[:,:,lev50].squeeze()
sortedobs = np.sort(obs50[:,1],axis=0)
sortedobsargs = np.argsort(obs50[:,1],axis=0)
strength = np.nanmean(obspchf[sortedobsargs[-5:],:,:] - obspchf[sortedobsargs[:5],:,:],axis=0)

### Loop through echam
strength_echam = np.empty((len(echam_obsu),months.shape[0],levels.shape[0]))
for i in range(len(strength_echam)):
    echam50 = echam_obsu[:,:,:,lev50]
    sorted_echam = np.sort(echam50[i,:,1],axis=0)
    sortedargs_echam = np.argsort(echam50[i,:,1],axis=0)
    strength_echam[i] = np.nanmean(echam_obsu[i,sortedargs_echam[-5:],:,:] - echam_obsu[i,sortedargs_echam[:5],:,:],axis=0)
strength_echamm = np.nanmean(strength_echam,axis=0)
diff_strength_echam = strength_echamm - strength

### Loop through cam
strength_cam = np.empty((len(cam_obsu),months.shape[0],levels.shape[0]))
for i in range(len(strength_cam)):
    cam50 = cam_obsu[:,:,:,lev50]
    sorted_cam = np.sort(cam50[i,:,1],axis=0)
    sortedargs_cam = np.argsort(cam50[i,:,1],axis=0)
    strength_cam[i] = np.nanmean(cam_obsu[i,sortedargs_cam[-5:],:,:] - cam_obsu[i,sortedargs_cam[:5],:,:],axis=0)
strength_camm = np.nanmean(strength_cam,axis=0)
diff_strength_cam = strength_camm - strength

### Loop through spear
strength_spear = np.empty((len(spear_obsu),months.shape[0],levels.shape[0]))
for i in range(len(strength_spear)):
    spear50 = spear_obsu[:,:,:,lev50]
    sorted_spear = np.sort(spear50[i,:,1],axis=0)
    sortedargs_spear = np.argsort(spear50[i,:,1],axis=0)
    strength_spear[i] = np.nanmean(spear_obsu[i,sortedargs_spear[-5:],:,:] - spear_obsu[i,sortedargs_spear[:5],:,:],axis=0)
strength_spearm = np.nanmean(strength_spear,axis=0)
diff_strength_spear = strength_spearm - strength

###############################################################################
###############################################################################
###############################################################################      
assembleForPlotting = [strength_echamm,strength,diff_strength_echam,
                        strength_camm,strength,diff_strength_cam,
                        strength_spearm,strength,diff_strength_spear]
titlelabels = ['AMIP-OBS','ERA5','DIFFERENCE']
labelmodels = np.repeat(['ECHAM5','CAM5','SPEAR'],3)

###########################################################################
###########################################################################
###########################################################################
##### Plot profiles
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
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
        
### Set limits for contours and colorbars6
limit = np.arange(-20,20.1,0.5)
barlim = np.arange(-20,21,5)
cmap = cmocean.cm.balance
label = r'\textbf{Composite of February-U50 of [U] m/s}'
zscale = np.array([1000,700,500,300,200,100,50,30,10])
        
fig = plt.figure(figsize=(10,4))

### Create plot
for p in range(len(assembleForPlotting)):
    ax1 = plt.subplot(3,3,p+1)
    ax1.spines['top'].set_color('dimgrey')
    ax1.spines['right'].set_color('dimgrey')
    ax1.spines['bottom'].set_color('dimgrey')
    ax1.spines['left'].set_color('dimgrey')
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['right'].set_linewidth(2)
    ax1.spines['top'].set_linewidth(2)
    ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='dimgrey')   
    ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='dimgrey')
    plt.gca().axes.get_xaxis().set_visible(True)
    plt.gca().axes.get_yaxis().set_visible(True)
    
    if any([p==0,p==3,p==6]):
        plt.ylabel(r'\textbf{Pressure [hPa]}',color='k',fontsize=7)
    
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    
    var = assembleForPlotting[p]
    
    ### Plot contours
    cs = plt.contourf(month2,lev2,var.transpose(),limit,extend='both')
    # cs1 = plt.contour(timeq,levq,climowinter.transpose(),np.arange(-10,35,2),
    #                   linewidths=2,colors='dimgrey')
    
    cs.set_cmap(cmap)
    
    plt.gca().invert_yaxis()
    plt.yscale('log')
    
    if any([p==6,p==7,p==8]):
        plt.xticks(np.arange(0,11+1,1),monthq,fontsize=4)
    else:
        plt.xticks([])
        
    if p < 3:
        plt.title(r'\textbf{%s}' % titlelabels[p],color='k',fontsize=12)
    
    if any([p==0,p==3,p==6]):
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    else:
        plt.yticks([])
        
    if any([p==2,p==5,p==8]):
        plt.text(11.6,100,r'\textbf{%s}' % labelmodels[p],fontsize=15,rotation=270,
                  color='dimgrey',ha='center',va='center')
    
    plt.xlim([0,11])
    plt.ylim([1000,10])
    plt.minorticks_off()

###########################################################################
plt.tight_layout()
cbar_ax = fig.add_axes([0.32,0.08,0.4,0.03])                    
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(label,fontsize=11,color='dimgrey',labelpad=1.4)  

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001,labelsize=7)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(bottom=0.17,hspace=0.12,wspace=0.08)    
plt.savefig(directoryfigure + 'HeatMap_FEB-U50-STRONGWEAK-PV_ZonalMeanZonalWindBias-AMIPS-obs_validationSPEAR.png',dpi=300)
print('Completed: Script done!')
