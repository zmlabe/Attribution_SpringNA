"""
Calculate Alaska Blocking Index

Reference : Ballinger et al. 2022 (IJOC)
Author    : Zachary M. Labe
Date      : 16 November 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
import calc_DetrendData as DT

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Create months to loop through
monthsall = ['January','February','March','April','May','JFM','FM','FMA','DJF','AMJ','JJA']

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

### Parameters
for i in range(len(monthsall)):
    variq = 'Z500'
    sliceperiod = monthsall[i]
    yearsall = np.arange(1979,2021+1,1)
    sliceshape = 3
    slicenan = 'nan'
    addclimo = True
    yearmin = 1979
    yearmax = 2020

    ### Read data
    latobs,lonobs,lev,varn = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                                 yearsall,sliceshape,addclimo,slicenan,'surface')
    lon2,lat2 = np.meshgrid(lonobs,latobs)
    
    ### Read only 1979-2020
    yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
    years = yearsall[yearq]
    var = varn[yearq,:,:]
    
    ### Calculate anomalies
    anoms = calc_anomalies(years,var)
    
    ### Detrend data
    # vardt = DT.detrendDataR(anoms,'surface','monthly')
    vardt = anoms
    
    ### Calculate ABI
    lonq = np.where((lonobs >=180) & (lonobs <=235))[0]
    latq = np.where((latobs >=55) & (latobs <=75))[0]
    anomlon = vardt[:,:,lonq]
    anoms = anomlon[:,latq,:]
    lat2sq = lat2[latq,:]
    lat2s = lat2sq[:,lonq]
    ABI = UT.calc_weightedAve(anoms,lat2s)
    
    ### Save index
    directoryoutput = '/work/Zachary.Labe/Data/ClimateIndices/ABI/'
    np.savetxt(directoryoutput + 'ABI_%s_%s_%s-%s.txt' % (variq,sliceperiod,
                                                                    yearmin,yearmax),
               np.c_[years,ABI])
    print('\n========Calculated Alaska Blocking Index=======\n')
