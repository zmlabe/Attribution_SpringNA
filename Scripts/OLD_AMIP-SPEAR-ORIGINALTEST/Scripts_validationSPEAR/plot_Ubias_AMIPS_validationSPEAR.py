"""
Plot bias in zonal-mean zonal wind of the AMIP experiments

Author    : Zachary M. Labe
Date      : 8 July 2022
"""

import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_FACTS_AMIPS as AM
import scipy.stats as sts
import itertools
import read_ERA5_monthly1x1 as ER

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
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
    elif data.ndim == 5:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:,:]
    elif data.ndim == 6:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:,:,:]
    
    return anoms

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

### Read in model data
meananom = []
diff = []
for s in range(len(scenario)):
    meananome = []
    diffanome = []
    for m in range(len(model[s])):
        lat1q,lon1q,lev1q,varq,yearsq = AM.read_FACTS_Experi(scenario[s],model[s][m],
                                                     variq,slicemonth,6,
                                                     slicenan,level)
        
        yearAllDataq = np.where((yearsq >= yearAllData.min()) & (yearsq <= yearAllData.max()))[0]
        varq = varq[:,yearAllDataq,:,:,:,:]
        print(yearAllDataq)
        print(lev1q)
        
        ### Calculate regional mean
        mean = calc_regionalAve('PCH',lat1q,lon1q,varq)
        
        ### Calculate bias
        bias = mean - obspchf[np.newaxis,:,:,:]
        
        meananome.append(mean)
        diffanome.append(bias)
    meananom.append(meananome)
    diff.append(diffanome)
    
### Save output
np.savez(directoryoutput + 'AMIP_SPEARval_MonthlyData_U-PCH.npz',u=meananom,levels=lev1q)
np.savez(directoryoutput + 'AMIP_SPEARval_MonthlyData_U-PCH_Bias.npz',bias=diff,levels=lev1q)
