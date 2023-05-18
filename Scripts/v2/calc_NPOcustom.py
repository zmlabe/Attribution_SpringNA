"""
Calculate North Pacific Oscillation for observations using custom index

Author    : Zachary M. Labe
Date      : 30 November 2022
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
    variq = 'SLP'
    sliceperiod = monthsall[i]
    yearsall = np.arange(1979,2021+1,1)
    sliceshape = 3
    slicenan = 'nan'
    addclimo = True
    yearmin = 1979
    yearmax = 2020

    ### Read data
    lat,lon,lev,varn = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                                 yearsall,sliceshape,addclimo,slicenan,'surface')
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read only 1979-2020
    yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
    years = yearsall[yearq]
    var = varn[yearq,:,:]
    
    ### Calculate anomalies
    anoms = calc_anomalies(years,var)
    
    ### Detrend data
    data = DT.detrendDataR(anoms,'surface','monthly')
    # data = anoms
    
    ### Calculate PNA
    ############################################################
    ### North SLP box
    latq_slp1 = np.where((lat >= 55) & (lat <= 70))[0]
    lonq_slp1 = np.where((lon >= 195) & (lon <= 230))[0]
    slp1_a = data[:,latq_slp1,:]
    slp1_b = slp1_a[:,:,lonq_slp1]
    lon2_slp1,lat2_slp1 = np.meshgrid(lon[lonq_slp1],lat[latq_slp1])
    slp_1final = UT.calc_weightedAve(slp1_b,lat2_slp1)  
    ############################################################
    ### South SLP box
    latq_slp2 = np.where((lat >= 10) & (lat <= 30))[0]
    lonq_slp2 = np.where((lon >= 183) & (lon <= 213))[0]
    slp2_a = data[:,latq_slp2,:]
    slp2_b = slp2_a[:,:,lonq_slp2]
    lon2_slp2,lat2_slp2 = np.meshgrid(lon[lonq_slp2],lat[latq_slp2])
    slp_2final = UT.calc_weightedAve(slp2_b,lat2_slp2)
        
    ############################################################
    ### Calculate NPOcustom
    NPOcustom_raw = slp_1final - slp_2final
    
    ### Save index
    directoryoutput = '/work/Zachary.Labe/Data/ClimateIndices/NPOcustom/'
    np.savetxt(directoryoutput + 'NPOcustom_%s_%s_%s-%s_detrended.txt' % (variq,sliceperiod,
                                                                    yearmin,yearmax),
               np.c_[years,NPOcustom_raw])
    print('\n========Calculated NPOcustom=======\n')
