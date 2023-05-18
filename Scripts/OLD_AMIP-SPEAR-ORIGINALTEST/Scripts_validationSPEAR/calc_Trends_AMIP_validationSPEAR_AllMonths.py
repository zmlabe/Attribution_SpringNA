"""
Plot histogram of trends over different time periods compared to SMILEs

Author    : Zachary M. Labe
Date      : 23 May 2022
"""

import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_FACTS_AMIPS as AM
import scipy.stats as sts
import itertools

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_obs_rf','spear']
model = ['ECHAM5','ESRL-CAM5','spear']
model_spear = ['spear']
variq = 'Z30'
slicemonth= 'none'
slicenan = 'nan'
addclimo = True
trendper = 'all'
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
        climold = np.nanmean(data[:,yearqold,:,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:,:]
    
    return anoms

def calc_regionalAve(regionn,lat1,lon1,data):
    if np.nanmax(lon1) < 200:
        print(ValueError('SOMETHING IS WRONG WITH THE LONGITUDE SCALE!'))
        sys.exit()
        
    if regionn == 'CANMID':
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
    elif regionn == 'ARCTIC':
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
        
    ### Calculate timeseries
    mean = UT.calc_weightedAve(meanbox,lat2q)

    return mean,lat1a,lat2q,lon1a,lon2q

def trendPeriods(years,trendper,data):
        yearModel = years
        if trendper == 'all':
            yearmin = 1979
            yearmax = yearModel.max()
        elif trendper == 'recent':
            yearmin = 2000
            yearmax = yearModel.max() 
        print(yearmin,yearmax)
        trend = UT.linearTrend(data,yearModel,'surface',yearmin,yearmax)*10
        
        return trend

def readData(scenario,model,slicemonth,yearAllData,trendper):
    lat1q,lon1q,varq,yearsq = AM.read_FACTS_Experi(scenario,model,
                                                 variq,slicemonth,5,
                                                 slicenan)
    
    ### Slice time period only from 1979-2019
    yearAllDataq = np.where((yearsq >= yearAllData.min()) & (yearsq <= yearAllData.max()))[0]
    varq = varq[:,yearAllDataq,:,:]
    
    ### Calculate anomalies
    anomsq = calc_anomalies(yearAllData,varq)
    
    ### Calculate regional mean
    mean,lat1r,lat2r,lon1r,lon2r = calc_regionalAve('ARCTIC',lat1q,lon1q,anomsq)
    
    ### Calculate trend for each month
    trend = []
    for i in range(anomsq.shape[2]):
        trendq = trendPeriods(yearAllData,trendper,anomsq[:,:,i,:,:])
        trend.append(trendq)
        print('Month is -------> %s!' % i)
        
    return mean,trend,lat1q,lon1q,lat1r,lon1r

### Read in data
mean_echam,trend_echam,lat1q_echam,lon1q_echam,lat1r_echam,lon1r_echam = readData(scenario[0],model[0],slicemonth,yearAllData,trendper)
mean_cam,trend_cam,lat1q_cam,lon1q_cam,lat1r_cam,lon1r_cam = readData(scenario[0],model[1],slicemonth,yearAllData,trendper)
mean_spear,trend_spear,lat1q_spear,lon1q_spear,lat1r_spear,lon1r_spear = readData(scenario[1],model[2],slicemonth,yearAllData,trendper)
          
### Save metadata and trends for AMIPS
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/'
np.savez(directoryoutput + 'AMIP_SPEARval_MonthlyData_%s_%s.npz' % ('AllMonths',variq),
         trend_echam = trend_echam,
         trend_cam = trend_cam,
         trend_spear = trend_spear,
         mean_echam = mean_echam,
         mean_cam = mean_cam,
         mean_spear = mean_spear,
         lat1_echam = lat1q_echam,
         lon1_echam = lon1q_echam,
         lat1_cam = lat1q_cam,
         lon1_cam = lon1q_cam,
         lat1_spear = lat1q_spear,
         lon1_spear = lon1q_spear,
         lat1r_echam = lat1r_echam,
         lon1r_echam = lon1r_echam,
         lat1r_cam = lat1r_cam,
         lon1r_cam = lon1r_cam,
         lat1r_spear = lat1r_spear,
         lon1r_spear = lon1r_spear)
