"""
Calculate North Pacific Oscillation for observations using custom index for 20CRv3

Author    : Zachary M. Labe
Date      : 30 November 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_20CRv3_monthly as CR
import calc_DetrendData as DT

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/20CRv3/' 

### Create months to loop through
monthsall = ['January','February','March','April','May','JFM','FM','FMA','AMJ','JJA']

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
    yearsall = np.arange(1921,2015+1,1)
    sliceshape = 3
    slicenan = 'nan'
    addclimo = True
    yearmin = 1921
    yearmax = 2015
    
    ### Read data
    lat,lon,var = CR.read_20CRv3_monthly(variq,directorydata,sliceperiod,
                                    yearsall,sliceshape,addclimo,slicenan)
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate anomalies
    anoms = np.asarray(calc_anomalies(yearsall,var))
    
    ### Detrend data
    data = DT.detrendDataR(anoms,'surface','monthly')
    # data = anoms
    
    ### Calculate PNA
    ############################################################
    ### North SLP box
    latq_slp1 = np.where((lat >= 55) & (lat <= 70))[0]
    lonq_slp1 = np.where((lon >= 175) & (lon <= 205))[0]
    slp1_a = data[:,latq_slp1,:]
    slp1_b = slp1_a[:,:,lonq_slp1]
    lon2_slp1,lat2_slp1 = np.meshgrid(lon[lonq_slp1],lat[latq_slp1])
    slp_1final = UT.calc_weightedAve(slp1_b,lat2_slp1)  
    ############################################################
    ### South SLP box
    latq_slp2 = np.where((lat >= 22) & (lat <= 35))[0]
    lonq_slp2 = np.where((lon >= 183) & (lon <= 213))[0]
    slp2_a = data[:,latq_slp2,:]
    slp2_b = slp2_a[:,:,lonq_slp2]
    lon2_slp2,lat2_slp2 = np.meshgrid(lon[lonq_slp2],lat[latq_slp2])
    slp_2final = UT.calc_weightedAve(slp2_b,lat2_slp2)
        
    ############################################################
    ### Calculate NPOcustomZ
    NPOcustomZ_raw = slp_1final - slp_2final
    
    ### Save index
    directoryoutput = '/work/Zachary.Labe/Data/ClimateIndices/NPOcustomZ/'
    np.savetxt(directoryoutput + 'NPOcustomZ-20CRv3_%s_%s_%s-%s_detrended.txt' % (variq,sliceperiod,
                                                                    yearmin,yearmax),
               np.c_[yearsall,NPOcustomZ_raw])
    print('\n========Calculated NPOcustomZ=======\n')
