"""
Plot time series of temperature over the region compared to large ensembles

Author    : Zachary M. Labe
Date      : 22 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthly as ERA
import read_SMILE_historical as SM
import read_LENS_historical as LE
import read_SPEAR_MED as SP
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ERA5/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'April'
years = np.arange(1979,2022+1,1)
slicenan = 'nan'
addclimo = True
newon = True
datareader = True

### Read data
if datareader == True:
    lat1,lon1,obs = ERA.read_ERA5_monthly(variq,'/work/Zachary.Labe/Data/ERA5/',sliceperiod,
                                    years,3,addclimo,
                                    slicenan,newon)
    lat1,lon1,gfdlc = SM.read_SMILEhistorical('/work/Zachary.Labe/Data/SMILE/','GFDL_CM3',variq,
                                                sliceperiod,4,
                                                slicenan,20)
    lat1,lon1,gfdlem = SM.read_SMILEhistorical('/work/Zachary.Labe/Data/SMILE/','GFDL_ESM2M',variq,
                                                sliceperiod,4,
                                                slicenan,20)
    lat1,lon1,lens = LE.read_LENShistorical('/work/Zachary.Labe/Data/LENS/monthly/',variq,
                                            sliceperiod,4,
                                            slicenan,40)
    lat1s,lon1s,spear = SP.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',variq,
                                            sliceperiod,4,
                                            slicenan,30,'all')
    
def calc_anomalies(years,data):
    """ 
    Calculate anomalies
    """
    
    ### Baseline - 1981-2010
    if data.ndim == 3:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[yearqold,:,:],axis=0)
        anoms = data - climold
    elif data.ndim == 4:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:]
    
    return anoms

def calc_regionalAve(regionn,lat1,lon1,data):
    if regionn == 'CANMID':
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
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
        
    ### Calculate timeseries
    mean = UT.calc_weightedAve(meanbox,lat2q)
    
    return mean,lat1a,lat2q,lon1a,lon2q

obs_anom = calc_anomalies(years,obs)
gfdlc_anom = calc_anomalies(years,gfdlc)
gfdlem_anom = calc_anomalies(years,gfdlem)
lens_anom = calc_anomalies(years,lens)
spear_anom = calc_anomalies(years,spear)

obs_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,obs_anom)
gfdlc_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,gfdlc_anom)
gfdlem_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,gfdlem_anom)
lens_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,lens_anom)
spear_reg,latrs1,latrs2,lonrs1,lonrs2 = calc_regionalAve('CANMID',lat1s,lon1s,spear_anom)

### Calculate statistics on observations
slope, intercept, r, p, se = sts.linregress(years[:-3],obs_reg[:-3]) # 1979-2019
trendline = slope*years[:-3] + intercept

slopeRE, interceptRE, rRE, pRE, seRE = sts.linregress(years[-20:],obs_reg[-20:]) # 2000-2019
trendlineRE = slopeRE*years[-20:] + interceptRE

### Decadal trend
dectrend = slope*10
dectrendRE = slopeRE*10

###############################################################################
###############################################################################
###############################################################################               
### Plot Figure
### Adjust axes in time series plots 
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
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',labelsize=6,pad=1.5)
ax.tick_params(axis='y',labelsize=6,pad=1.5)

plt.plot(years[:-3],obs_reg[:-3],linewidth=2.3,color='maroon',alpha=1,clip_on=False,
         marker='o',markersize=6,label='ERA5',zorder=6)
for i in range(len(lens_reg)):
    plt.plot(years[:-3],lens_reg[i,:-3],color='darkgrey',alpha=0.5,clip_on=False,
             linewidth=0.7,label='CESM1-LE')
    if i == 0:
        leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
              bbox_to_anchor=(0.5,0.06),fancybox=True,ncol=4,frameon=False,
              handlelength=1,handletextpad=0.5)
plt.plot(years[:-3],np.nanmean(lens_reg[:,:-3],axis=0),color='dimgrey',alpha=1,
         linestyle='--',dashes=(1,0.3),linewidth=2.3)

plt.plot(years[:-3],trendline,linewidth=1,color='teal',clip_on=False)
plt.plot(years[-23:-3],trendlineRE,linewidth=3,color='darkblue',clip_on=False)

plt.xticks(np.arange(1980,2040,10),np.arange(1980,2040,10))
plt.yticks(np.round(np.arange(-8,8.1,1),2),np.round(np.arange(-8,8.1,1),2))
plt.xlim([1979,int(years[-1])-2])
plt.ylim([-7,7])

plt.ylabel(r'\textbf{$\bf{^\circ}$C}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{APRIL TEMPERATURE ANOMALIES}',
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'April_meanT_1979-2019_CESM1.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################               
### Plot Figure
### Adjust axes in time series plots 
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
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',labelsize=6,pad=1.5)
ax.tick_params(axis='y',labelsize=6,pad=1.5)

plt.plot(years[:-3],obs_reg[:-3],linewidth=2.3,color='maroon',alpha=1,clip_on=False,
         marker='o',markersize=6,label='ERA5',zorder=6)
for i in range(len(gfdlc_reg)):
    plt.plot(years[:-3],gfdlc_reg[i,:-3],color='darkgrey',alpha=0.5,clip_on=False,
             linewidth=0.7,label='GFDL-CM3')
    if i == 0:
        leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
              bbox_to_anchor=(0.5,0.06),fancybox=True,ncol=4,frameon=False,
              handlelength=1,handletextpad=0.5)
plt.plot(years[:-3],np.nanmean(gfdlc_reg[:,:-3],axis=0),color='dimgrey',alpha=1,
         linestyle='--',dashes=(1,0.3),linewidth=2.3)
plt.plot(years[:-3],trendline,linewidth=1,color='teal',clip_on=False)
plt.plot(years[-23:-3],trendlineRE,linewidth=3,color='darkblue',clip_on=False)

plt.xticks(np.arange(1980,2040,10),np.arange(1980,2040,10))
plt.yticks(np.round(np.arange(-8,8.1,1),2),np.round(np.arange(-8,8.1,1),2))
plt.xlim([1979,int(years[-1])-2])
plt.ylim([-7,7])

plt.ylabel(r'\textbf{$\bf{^\circ}$C}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{APRIL TEMPERATURE ANOMALIES}',
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'April_meanT_1979-2019_GFDLCM3.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################               
### Plot Figure
### Adjust axes in time series plots 
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
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',labelsize=6,pad=1.5)
ax.tick_params(axis='y',labelsize=6,pad=1.5)

plt.plot(years[:-3],obs_reg[:-3],linewidth=2.3,color='maroon',alpha=1,clip_on=False,
         marker='o',markersize=6,label='ERA5',zorder=6)
for i in range(len(gfdlem_reg)):
    plt.plot(years[:-3],gfdlem_reg[i,:-3],color='darkgrey',alpha=0.5,clip_on=False,
             linewidth=0.7,label='GFDL-ESM2M')
    if i == 0:
        leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
              bbox_to_anchor=(0.5,0.06),fancybox=True,ncol=4,frameon=False,
              handlelength=1,handletextpad=0.5)
plt.plot(years[:-3],np.nanmean(gfdlem_reg[:,:-3],axis=0),color='dimgrey',alpha=1,
         linestyle='--',dashes=(1,0.3),linewidth=2.3)        
        
plt.plot(years[:-3],trendline,linewidth=1,color='teal',clip_on=False)
plt.plot(years[-23:-3],trendlineRE,linewidth=3,color='darkblue',clip_on=False)

plt.xticks(np.arange(1980,2040,10),np.arange(1980,2040,10))
plt.yticks(np.round(np.arange(-8,8.1,1),2),np.round(np.arange(-8,8.1,1),2))
plt.xlim([1979,int(years[-1])-2])
plt.ylim([-7,7])

plt.ylabel(r'\textbf{$\bf{^\circ}$C}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{APRIL TEMPERATURE ANOMALIES}',
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'April_meanT_1979-2019_GFDLEM.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################               
### Plot Figure
### Adjust axes in time series plots 
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
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',labelsize=6,pad=1.5)
ax.tick_params(axis='y',labelsize=6,pad=1.5)

plt.plot(years[:-3],obs_reg[:-3],linewidth=2.3,color='maroon',alpha=1,clip_on=False,
         marker='o',markersize=6,label='ERA5',zorder=6)
for i in range(len(spear_reg)):
    plt.plot(years[:-3],spear_reg[i,:-3],color='darkgrey',alpha=0.5,clip_on=False,
             linewidth=0.7,label='SPEAR_MED')
    if i == 0:
        leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
              bbox_to_anchor=(0.5,0.06),fancybox=True,ncol=4,frameon=False,
              handlelength=1,handletextpad=0.5)
plt.plot(years[:-3],np.nanmean(spear_reg[:,:-3],axis=0),color='dimgrey',alpha=1,
         linestyle='--',dashes=(1,0.3),linewidth=2.3)        
        
plt.plot(years[:-3],trendline,linewidth=1,color='teal',clip_on=False)
plt.plot(years[-23:-3],trendlineRE,linewidth=3,color='darkblue',clip_on=False)

plt.xticks(np.arange(1980,2040,10),np.arange(1980,2040,10))
plt.yticks(np.round(np.arange(-8,8.1,1),2),np.round(np.arange(-8,8.1,1),2))
plt.xlim([1979,int(years[-1])-2])
plt.ylim([-7,7])

plt.ylabel(r'\textbf{$\bf{^\circ}$C}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{APRIL TEMPERATURE ANOMALIES}',
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'April_meanT_1979-2019_SPEAR.png',dpi=300)
