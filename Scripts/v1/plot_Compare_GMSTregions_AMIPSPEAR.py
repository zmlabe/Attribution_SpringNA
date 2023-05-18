"""
Compare regional mean temperatures for different SPEAR experiments
 
Author    : Zachary M. Labe
Date      : 2 September 2022
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
import read_SPEAR_MED_Scenario as sc
import read_FACTS_AMIPS as FA
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/AMIPS/GMSTregions/' 

### Parameters
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
slicemonth = 'April'
slicenan = 'nan'
years = np.arange(1979,2020+1,1)

###############################################################################

def remove_ocean(data):
    """
    Masks out the ocean for land_only == True
    """
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Read in land mask
    directorydata = '/work/Zachary.Labe/Data/masks/'
    filename = 'land_maskcoarse_SPEAR_MED.nc'
    datafile = Dataset(directorydata + filename)
    mask = np.asarray(datafile.variables['land_mask'][:])
    datafile.close()
    
    ### Mask out model and observations
    datamask = data * mask
    
    ### Check for floats
    datamask[np.where(datamask==0.)] = np.nan
    
    return datamask

###############################################################################

def remove_land(data):
    """
    Masks out the ocean for ocean_only == True
    """
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Read in ocean mask
    directorydata = '/work/Zachary.Labe/Data/masks/'
    filename = 'ocean_maskcoarse_SPEAR_MED.nc'
    datafile = Dataset(directorydata + filename)
    mask = np.asarray(datafile.variables['ocean_mask'][:])
    datafile.close()
    
    ### Mask out model and observations
    datamask = data * mask
    
    ### Check for floats
    datamask[np.where(datamask==0.)] = np.nan
    
    return datamask

def calcAverages(data,lat,lon,region,maskocean):
    ### Calculate anomalies
    if data.ndim == 3:
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        clim = np.nanmean(data[yearq,:,:],axis=0)
        data = data - clim
    elif data.ndim == 4:
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        clim = np.nanmean(data[:,yearq,:,:],axis=1)
        data = data - clim[:,np.newaxis,:,:]
    else:
        print(ValueError('WRONG DIMENSIONS OF DATA!'))
        sys.exit()
    
    ### Only plot land
    if maskocean == True:
        data = remove_ocean(data)
    
    ### Calculate averages
    if data.ndim == 3:
        if region == 'Globe':
            latnew = lat
            lon2,lat2 = np.meshgrid(lon,latnew)
            datanew = data[:,:,:]
            ave = UT.calc_weightedAve(datanew,lat2)
        elif region == 'NH':
            latq = np.where((lat > 0))[0]
            latnew = lat[latq]
            lon2,lat2 = np.meshgrid(lon,latnew)
            datanew = data[:,latq,:]
            ave = UT.calc_weightedAve(datanew,lat2)
        elif region == 'SH':
            latq = np.where((lat < 0))[0]
            latnew = lat[latq]
            lon2,lat2 = np.meshgrid(lon,latnew)
            datanew = data[:,latq,:]
            ave = UT.calc_weightedAve(datanew,lat2)
        else:
            print(ValueError('WRONG REGION SELECTED!'))
            sys.exit()
    elif data.ndim == 4:
        if region == 'Globe':
            latnew = lat
            lon2,lat2 = np.meshgrid(lon,latnew)
            datanew = data[:,:,:,:]
            ave = UT.calc_weightedAve(datanew,lat2)
        elif region == 'NH':
            latq = np.where((lat > 0))[0]
            latnew = lat[latq]
            lon2,lat2 = np.meshgrid(lon,latnew)
            datanew = data[:,:,latq,:]
            ave = UT.calc_weightedAve(datanew,lat2)
        elif region == 'SH':
            latq = np.where((lat < 0))[0]
            latnew = lat[latq]
            lon2,lat2 = np.meshgrid(lon,latnew)
            datanew = data[:,:,latq,:]
            ave = UT.calc_weightedAve(datanew,lat2)
        else:
            print(ValueError('WRONG REGION SELECTED!'))
            sys.exit()
    else:
        print(ValueError('WRONG DIMENSIONS OF DATA!'))
        sys.exit()
        
    return ave


### Read in Observations
latobs,lonobs,levobs,var = ER.read_ERA5_monthlyMEDS(variq,'/work/Zachary.Labe/Data/',slicemonth,years,3,True,slicenan,'surface')
era = var[:-1,:,:] # through 2020

### Read in AMIP-SPEAR
lat,lon,lev,amipq,yearsf = FA.read_FACTS_Experi('spear_obs_rf','SPEAR',variq,slicemonth,4,slicenan,'surface')
yearamip = np.arange(1901,2020+1,1)
yearamipq = np.where((yearamip >= 1979) & (yearamip <= 2020))[0]
amip = amipq[:,yearamipq,:,:]

### Read in SPEAR_MED
lat1s,lon1s,spear = SP.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',variq,slicemonth,4,slicenan,30,'all')
spearm = np.nanmean(spear,axis=0)

### Read in SPEAR_MED-SSP245
lat1s,lon1s,spearsp = sc.read_SPEAR_MED_Scenario('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/','SSP245',variq,slicemonth,4,slicenan,30,'all')
spearspm = np.nanmean(spearsp,axis=0)

### Calculate averages
era_globe = calcAverages(era,lat,lon,'Globe',False)
amip_globe = calcAverages(amip,lat,lon,'Globe',False)
spear_globe = calcAverages(spearm,lat,lon,'Globe',False)
spearsp_globe = calcAverages(spearspm,lat,lon,'Globe',False)

era_globeland = calcAverages(era,lat,lon,'Globe',True)
amip_globeland = calcAverages(amip,lat,lon,'Globe',True)
spear_globeland = calcAverages(spearm,lat,lon,'Globe',True)
spearsp_globeland = calcAverages(spearspm,lat,lon,'Globe',True)

era_NH = calcAverages(era,lat,lon,'NH',False)
amip_NH = calcAverages(amip,lat,lon,'NH',False)
spear_NH = calcAverages(spearm,lat,lon,'NH',False)
spearsp_NH = calcAverages(spearspm,lat,lon,'NH',False)

era_NHland = calcAverages(era,lat,lon,'NH',True)
amip_NHland = calcAverages(amip,lat,lon,'NH',True)
spear_NHland = calcAverages(spearm,lat,lon,'NH',True)
spearsp_NHland = calcAverages(spearspm,lat,lon,'NH',True)

era_SH = calcAverages(era,lat,lon,'SH',False)
amip_SH = calcAverages(amip,lat,lon,'SH',False)
spear_SH = calcAverages(spearm,lat,lon,'SH',False)
spearsp_SH = calcAverages(spearspm,lat,lon,'SH',False)

era_SHland = calcAverages(era,lat,lon,'SH',True)
amip_SHland = calcAverages(amip,lat,lon,'SH',True)
spear_SHland = calcAverages(spearm,lat,lon,'SH',True)
spearsp_SHland = calcAverages(spearspm,lat,lon,'SH',True)

slope_era_globe,intercept,r,p,se = sts.linregress(years,era_globe)
slope_spear_globe,intercept,r,p,se = sts.linregress(years,spear_globe)
slope_spearsp_globe,intercept,r,p,se = sts.linregress(years,spearsp_globe)
slope_amip_globe,intercept,r,p,se = sts.linregress(years,np.nanmean(amip_globe,axis=0))

slope_era_globeland,intercept,r,p,se = sts.linregress(years,era_globeland)
slope_spear_globeland,intercept,r,p,se = sts.linregress(years,spear_globeland)
slope_spearsp_globeland,intercept,r,p,se = sts.linregress(years,spearsp_globeland)
slope_amip_globeland,intercept,r,p,se = sts.linregress(years,np.nanmean(amip_globeland,axis=0))

slope_era_NHland,intercept,r,p,se = sts.linregress(years,era_NHland)
slope_spear_NHland,intercept,r,p,se = sts.linregress(years,spear_NHland)
slope_spearsp_NHland,intercept,r,p,se = sts.linregress(years,spearsp_NHland)
slope_amip_NHland,intercept,r,p,se = sts.linregress(years,np.nanmean(amip_NHland,axis=0))

slope_era_SHland,intercept,r,p,se = sts.linregress(years,era_SHland)
slope_spear_SHland,intercept,r,p,se = sts.linregress(years,spear_SHland)
slope_spearsp_SHland,intercept,r,p,se = sts.linregress(years,spearsp_SHland)
slope_amip_SHland,intercept,r,p,se = sts.linregress(years,np.nanmean(amip_SHland,axis=0))

###########################################################################
###########################################################################
###########################################################################
fig = plt.figure(figsize=(8,4))
ax = plt.subplot(121)

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

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')
ax.tick_params(axis='x',labelsize=7,pad=4)
ax.tick_params(axis='y',labelsize=7,pad=4)

plt.plot(years,amip_globe.transpose(),linestyle='-',linewidth=0.3,color='deepskyblue',
         alpha=0.4,clip_on=False)
plt.plot(years,era_globe,linestyle='--',linewidth=1.5,color='k',
         dashes=(1,0.3),clip_on=False,label=r'\textbf{ERA5}')
plt.plot(years,spearsp_globe,linestyle='-',linewidth=2,color='lightcoral',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP245}')
plt.plot(years,spear_globe,linestyle='-',linewidth=2,color='crimson',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP585}')
plt.plot(years,np.nanmean(amip_globe,axis=0),linestyle='-',linewidth=1.5,color='darkblue',
         alpha=1,clip_on=False,label=r'\textbf{SPEAR-MED-AMIP}')

leg = plt.legend(shadow=False,fontsize=9,loc='upper center',
                  bbox_to_anchor=(0.5,1.02),fancybox=True,ncol=2,frameon=False,
                  handlelength=1,handletextpad=0.3)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
    
plt.text(2008,-0.5,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_era_globe*10,2)),
         color='k',fontsize=10)
plt.text(2008,-0.6,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spearsp_globe*10,2)),
         color='lightcoral',fontsize=10)
plt.text(2008,-0.7,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spear_globe*10,2)),
         color='crimson',fontsize=10)
plt.text(2008,-0.8,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_amip_globe*10,2)),
         color='darkblue',fontsize=10)

plt.title(r'\textbf{GMST [land+ocean]}',color='dimgrey',fontsize=15)
plt.ylabel(r'\textbf{%s Temperature Anomaly [$^{\circ}$C] - 1981-2010}' % slicemonth,fontsize=7,color='dimgrey')
plt.yticks(np.arange(-20,17,0.2),map(str,np.round(np.arange(-20,17,0.2),2)))
plt.xticks(np.arange(1980,2020+1,10),map(str,np.arange(1980,2020+1,10)))
plt.xlim([1979,2020])   
plt.ylim([-0.8,1.2])

###########################################################################
ax = plt.subplot(122)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')
ax.tick_params(axis='x',labelsize=7,pad=4)
ax.tick_params(axis='y',labelsize=7,pad=4)

plt.plot(years,amip_globeland.transpose(),linestyle='-',linewidth=0.3,color='deepskyblue',
         alpha=0.4,clip_on=False)
plt.plot(years,era_globeland,linestyle='--',linewidth=1.5,color='k',
         dashes=(1,0.3),clip_on=False,label=r'\textbf{ERA5}')
plt.plot(years,spearsp_globeland,linestyle='-',linewidth=2,color='lightcoral',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP245}')
plt.plot(years,spear_globeland,linestyle='-',linewidth=2,color='crimson',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP585}')
plt.plot(years,np.nanmean(amip_globeland,axis=0),linestyle='-',linewidth=1.5,color='darkblue',
         alpha=1,clip_on=False,label=r'\textbf{SPEAR-MED-AMIP}')

leg = plt.legend(shadow=False,fontsize=9,loc='upper center',
                  bbox_to_anchor=(0.5,1.02),fancybox=True,ncol=2,frameon=False,
                  handlelength=1,handletextpad=0.3)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
    
plt.text(2008,-0.5,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_era_globeland*10,2)),
         color='k',fontsize=10)
plt.text(2008,-0.6,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spearsp_globeland*10,2)),
         color='lightcoral',fontsize=10)
plt.text(2008,-0.7,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spear_globeland*10,2)),
         color='crimson',fontsize=10)
plt.text(2008,-0.8,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_amip_globeland*10,2)),
         color='darkblue',fontsize=10)

plt.title(r'\textbf{GMST [land]}',color='dimgrey',fontsize=15)
plt.yticks(np.arange(-20,17,0.2),map(str,np.round(np.arange(-20,17,0.2),2)))
plt.xticks(np.arange(1980,2020+1,10),map(str,np.arange(1980,2020+1,10)))
plt.xlim([1979,2020])   
plt.ylim([-0.8,1.2])

plt.tight_layout()
plt.savefig(directoryfigure + 'GMST_SPEARAMIP_%s.png' % slicemonth,dpi=300)

###########################################################################
###########################################################################
###########################################################################
fig = plt.figure(figsize=(8,4))
ax = plt.subplot(121)

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

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')
ax.tick_params(axis='x',labelsize=7,pad=4)
ax.tick_params(axis='y',labelsize=7,pad=4)

plt.plot(years,amip_NHland.transpose(),linestyle='-',linewidth=0.3,color='deepskyblue',
         alpha=0.4,clip_on=False)
plt.plot(years,era_NHland,linestyle='--',linewidth=1.5,color='k',
         dashes=(1,0.3),clip_on=False,label=r'\textbf{ERA5}')
plt.plot(years,spearsp_NHland,linestyle='-',linewidth=2,color='lightcoral',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP245}')
plt.plot(years,spear_NHland,linestyle='-',linewidth=2,color='crimson',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP585}')
plt.plot(years,np.nanmean(amip_NHland,axis=0),linestyle='-',linewidth=1.5,color='darkblue',
         alpha=1,clip_on=False,label=r'\textbf{SPEAR-MED-AMIP}')

leg = plt.legend(shadow=False,fontsize=9,loc='upper center',
                  bbox_to_anchor=(0.5,1.02),fancybox=True,ncol=2,frameon=False,
                  handlelength=1,handletextpad=0.3)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
    
plt.text(2008,-0.5,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_era_NHland*10,2)),
         color='k',fontsize=10)
plt.text(2008,-0.6,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spearsp_NHland*10,2)),
         color='lightcoral',fontsize=10)
plt.text(2008,-0.7,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spear_NHland*10,2)),
         color='crimson',fontsize=10)
plt.text(2008,-0.8,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_amip_NHland*10,2)),
         color='darkblue',fontsize=10)

plt.title(r'\textbf{NH [land]}',color='dimgrey',fontsize=15)
plt.ylabel(r'\textbf{%s Temperature Anomaly [$^{\circ}$C] - 1981-2010}' % slicemonth,fontsize=7,color='dimgrey')
plt.yticks(np.arange(-20,17,0.2),map(str,np.round(np.arange(-20,17,0.2),2)))
plt.xticks(np.arange(1980,2020+1,10),map(str,np.arange(1980,2020+1,10)))
plt.xlim([1979,2020])   
plt.ylim([-0.8,1.6])

###########################################################################
ax = plt.subplot(122)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')
ax.tick_params(axis='x',labelsize=7,pad=4)
ax.tick_params(axis='y',labelsize=7,pad=4)

plt.plot(years,amip_SHland.transpose(),linestyle='-',linewidth=0.3,color='deepskyblue',
         alpha=0.4,clip_on=False)
plt.plot(years,era_SHland,linestyle='--',linewidth=1.5,color='k',
         dashes=(1,0.3),clip_on=False,label=r'\textbf{ERA5}')
plt.plot(years,spearsp_SHland,linestyle='-',linewidth=2,color='lightcoral',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP245}')
plt.plot(years,spear_SHland,linestyle='-',linewidth=2,color='crimson',
         clip_on=False,label=r'\textbf{SPEAR-MED-SSP585}')
plt.plot(years,np.nanmean(amip_SHland,axis=0),linestyle='-',linewidth=1.5,color='darkblue',
         alpha=1,clip_on=False,label=r'\textbf{SPEAR-MED-AMIP}')

leg = plt.legend(shadow=False,fontsize=9,loc='upper center',
                  bbox_to_anchor=(0.5,1.02),fancybox=True,ncol=2,frameon=False,
                  handlelength=1,handletextpad=0.3)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
    
plt.text(2008,-0.5,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_era_SHland*10,2)),
         color='k',fontsize=10)
plt.text(2008,-0.6,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spearsp_SHland*10,2)),
         color='lightcoral',fontsize=10)
plt.text(2008,-0.7,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_spear_SHland*10,2)),
         color='crimson',fontsize=10)
plt.text(2008,-0.8,r'\textbf{%s$^{\circ}$/decade}' % (np.round(slope_amip_SHland*10,2)),
         color='darkblue',fontsize=10)

plt.title(r'\textbf{SH [land]}',color='dimgrey',fontsize=15)
plt.yticks(np.arange(-20,17,0.2),map(str,np.round(np.arange(-20,17,0.2),2)))
plt.xticks(np.arange(1980,2020+1,10),map(str,np.arange(1980,2020+1,10)))
plt.xlim([1979,2020])   
plt.ylim([-0.8,1.6])

plt.tight_layout()
plt.savefig(directoryfigure + 'HemisphereGMST_SPEARAMIP_%s.png' % slicemonth,dpi=300)
        
        
