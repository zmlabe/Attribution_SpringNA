"""
Calculate teleconnections for SPEAR-MED-NOAER

Author    : Zachary M. Labe
Date      : 21 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_SPEAR_MED_NOAER as SP
import calc_DetrendData as DTT

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/SND_BoxNA_AllModels/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','May','JFM','FM','FMA','AMJ','JJA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    detrend = True
    years = np.arange(1979,2019+1,1)
    
    indices = ['SND_BoxNA']
    indicesvari = ['SNOW']
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Functions for this script
    def calc_SND_BoxNA(data,lat,lon,years):
        """ Calculate Box of Interest over NA for any variables
        """
        ############################################################
        ### Calculated SND_BoxNA
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        latq = np.where((lat >= la1) & (lat <= la2))[0]
        lonq = np.where((lon >= lo1) & (lon <= lo2))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        SND_BoxNA_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize SND_BoxNA by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(SND_BoxNA_raw[yearq],axis=0)
        mean = np.mean(SND_BoxNA_raw[yearq],axis=0)
        SND_BoxNA = (SND_BoxNA_raw - mean)/std
        
        print('Calculated -----> SND_BoxNA Index!')
        return SND_BoxNA
    ###############################################################################
    ###############################################################################
    ###############################################################################
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
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE DIMENSIONS OF ANOMALIES!'))
            sys.exit()
        
        print('Finished calculating anomalies!')
        return anoms
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calcTeleconnections(indicesvari,indices,monq,years):
        spear_lat,spear_lon,spear = SP.read_SPEAR_MED_NOAER('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED_NOAER/monthly/',
                                              indicesvari,monq,4,'nan',12,'all')
        
    
        ### Calculate anomalies 
        spear_anom = calc_anomalies(years,spear) 
    
        ### Detrend data
        if detrend == True:
            spear_dt = DTT.detrendData(spear_anom,'surface','yearly')
        else:
            spear_dt = spear_anom
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'SND_BoxNA':
            ### Calculate indices for each ensemble member of SPEAR
            SND_BoxNA_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                SND_BoxNA_spear[ens,:] = calc_SND_BoxNA(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = SND_BoxNA_spear
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE INDEX SELECTED!'))
            sys.exit()
        return indexz_spear
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Begin loop for indices for calculations
    data_spear = np.empty((len(indices),12,years.shape[0]))
    for iii in range(len(indices)):
        data_spear[iii,:,:] = calcTeleconnections(indicesvari[iii],indices[iii],monq,years)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Save output 
    np.savez(directorydata + 'ClimateIndices_coupledSPEAR-NOAER-SND_BoxNA_%s_1979-2019.npz' % monq,
                 spear=data_spear,indices=indices,indicesvari=indicesvari,years=years)
