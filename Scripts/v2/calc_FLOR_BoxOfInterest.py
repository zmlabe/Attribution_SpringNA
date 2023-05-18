"""
Calculate teleconnections for FLOR

Author    : Zachary M. Labe
Date      : 5 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_FLOR as FL
import calc_DetrendData as DTT

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/FLOR/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','JFM','FM','FMA','AMJ','JJA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    detrend = True
    years = np.arange(1979,2020+1,1)
    
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
    
        ############################################################
        ### Normalize T2M_BoxNA by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(T2M_BoxNA_raw[yearq],axis=0)
        mean = np.mean(T2M_BoxNA_raw[yearq],axis=0)
        T2M_BoxNA = (T2M_BoxNA_raw - mean)/std
        
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
        FLOR_lat,FLOR_lon,FLOR = FL.read_FLOR('/work/Zachary.Labe/Data/',
                                              indicesvari,monq,4,'nan',
                                              30,'satellite')
        
    
        ### Calculate anomalies 
        FLOR_anom = calc_anomalies(years,FLOR) 
    
        ### Detrend data
        if detrend == True:
            FLOR_dt = DTT.detrendData(FLOR_anom,'surface','yearly')
        else:
            FLOR_dt = FLOR_anom
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'T2M_BoxNA':
            ### Calculate indices for each ensemble member of FLOR
            T2M_BoxNA_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                T2M_BoxNA_FLOR[ens,:] = calc_T2M_BoxNA(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = T2M_BoxNA_FLOR
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE INDEX SELECTED!'))
            sys.exit()
        return indexz_FLOR
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Begin loop for indices for calculations
    data_FLOR = np.empty((len(indices),30,years.shape[0]))
    for iii in range(len(indices)):
        data_FLOR[iii,:,:] = calcTeleconnections(indicesvari[iii],indices[iii],monq,years)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Save output 
    np.savez(directorydata + 'ClimateIndices_coupledFLOR-T2M_BoxNA_%s_1979-2020.npz' % monq,
                 FLOR=data_FLOR,indices=indices,indicesvari=indicesvari,years=years)
