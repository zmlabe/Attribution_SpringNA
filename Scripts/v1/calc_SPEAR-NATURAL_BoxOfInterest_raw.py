"""
Calculate teleconnections for SPEAR-MED-NATURAL and save raw anomalies

Author    : Zachary M. Labe
Date      : 15 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_SPEAR_MED_NATURAL as SPNAT

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/T2M_BoxNA_AllModels/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','JFM','FM','FMA','AMJ','JJA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    detrend = False
    years = np.arange(1979,2019+1,1)
    
    indices = ['T2M_BoxNA']
    indicesvari = ['T2M']
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Functions for this script
    def calc_T2M_BoxNA(data,lat,lon,years):
        """ Calculate Box of Interest over NA for any variables
        """
        ############################################################
        ### Calculated T2M_BoxNA
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        latq = np.where((lat >= la1) & (lat <= la2))[0]
        lonq = np.where((lon >= lo1) & (lon <= lo2))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        T2M_BoxNA_raw = UT.calc_weightedAve(var_b,lat2)  
        T2M_BoxNA = T2M_BoxNA_raw 
        
        print('Calculated -----> T2M_BoxNA Index!')
        return T2M_BoxNA
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
        spear_lat,spear_lon,spear = SPNAT.read_SPEAR_MED_NATURAL('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED_NATURAL/monthly/',
                                                                indicesvari,monq,4,
                                                                 'nan',30,
                                                                 'satellite')
        
        ### Calculate anomalies 
        spear_dt = calc_anomalies(years,spear) 
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'T2M_BoxNA':
            ### Calculate indices for each ensemble member of SPEAR
            T2M_BoxNA_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                T2M_BoxNA_spear[ens,:] = calc_T2M_BoxNA(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = T2M_BoxNA_spear
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE INDEX SELECTED!'))
            sys.exit()
        return indexz_spear
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Begin loop for indices for calculations
    data_spear = np.empty((len(indices),30,years.shape[0]))
    for iii in range(len(indices)):
        data_spear[iii,:,:] = calcTeleconnections(indicesvari[iii],indices[iii],monq,years)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Save output 
    np.savez(directorydata + 'ClimateIndices_coupledSPEAR-NATURAL-T2M_BoxNA_%s_1979-2019_raw.npz' % monq,
                 spear=data_spear,indices=indices,indicesvari=indicesvari,years=years)
