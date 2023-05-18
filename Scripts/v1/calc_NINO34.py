"""
Calculate NINO34

Author    : Zachary M. Labe
Date      : 29 June 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/' 

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
variq = 'SST'
sliceperiod = 'JFM'
yearsall = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
yearmin = 1979
yearmax = 2019

### Read data
latobs,lonobs,lev,varn = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                             yearsall,sliceshape,addclimo,slicenan,'surface')
lon2,lat2 = np.meshgrid(lonobs,latobs)

### Read only 1979-2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]
var = varn[yearq,:,:]

### Calculate anomalies
anoms = np.asarray(calc_anomalies(years,var))

### Detrend data
vardt = DT.detrendDataR(anoms,'surface','monthly')
# vardt = anoms

### Calculate SHI
lonq = np.where((lonobs >=190) & (lonobs <=240))[0]
latq = np.where((latobs >=-5) & (latobs <=5))[0]
anomlon = vardt[:,:,lonq]
anoms = anomlon[:,latq,:]
lat2sq = lat2[latq,:]
lat2s = lat2sq[:,lonq]
shi = UT.calc_weightedAve(anoms,lat2s)

### Save index
directoryoutput = '/work/Zachary.Labe/Data/ClimateIndices/NINO34/'
np.savetxt(directoryoutput + 'NINO34_%s_%s_%s-%s_detrended.txt' % (variq,sliceperiod,
                                                                yearmin,yearmax),
           np.c_[years,shi])
print('\n========Calculated NINO34 SST=======\n')
