"""
Calculate tropical teleconnections for SPEAR-MED-NOAER and do not detrend

Author    : Zachary M. Labe
Date      : 13 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_SPEAR_MED_NOAER as SP

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/coupledSPEAR-NOAER/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','May','JFM','FM','FMA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    detrend = False
    years = np.arange(1979,2019+1,1)
    
    indices = ['NINO34','WPPRECT','NECPSST']
    indicesvari = ['SST','PRECT','SST']
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Functions for this script
    def calc_WPPRECT(data,lat,lon,years):
        """ Calculate WPPRECT for any variables
        """
        ############################################################
        ### Calculated WPPRECT Index
        lonq = np.where((lon >=110) & (lon <=140))[0]
        latq = np.where((lat >=-5) & (lat <=25))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        WPPRECT_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize WPPRECT by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(WPPRECT_raw[yearq],axis=0)
        mean = np.mean(WPPRECT_raw[yearq],axis=0)
        WPPRECT = (WPPRECT_raw - mean)/std
        
        print('Calculated -----> WPPRECT Index!')
        return WPPRECT
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_NECPSST(data,lat,lon,years):
        """ Calculate NECPTSST for any variables
        """
        ############################################################
        ### Calculated NECPTSST Index
        lonq = np.where((lon >=180) & (lon <=245))[0]
        latq = np.where((lat >=0) & (lat <=25))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        NECPTSST_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize NECPTSST by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(NECPTSST_raw[yearq],axis=0)
        mean = np.mean(NECPTSST_raw[yearq],axis=0)
        NECPTSST = (NECPTSST_raw - mean)/std
        
        print('Calculated -----> NECPTSST Index!')
        return NECPTSST
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_NINO34(data,lat,lon,years):
        """ Calculate Nino 3.4 for any variables
        """
        ############################################################
        ### Calculated Nino 3.4 Index
        lonq = np.where((lon >=190) & (lon <=240))[0]
        latq = np.where((lat >=-5) & (lat <=5))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        NINO34_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize NINO34 by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(NINO34_raw[yearq],axis=0)
        mean = np.mean(NINO34_raw[yearq],axis=0)
        NINO34 = (NINO34_raw - mean)/std
        
        print('Calculated -----> NINO34 Index!')
        return NINO34
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
        spear_dt = calc_anomalies(years,spear) 
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'NINO34':
            ### Calculate indices for each ensemble member of SPEAR
            NINO34_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                NINO34_spear[ens,:] = calc_NINO34(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = NINO34_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'WPPRECT':
            ### Calculate indices for each ensemble member of SPEAR
            WPPRECT_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                WPPRECT_spear[ens,:] = calc_WPPRECT(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = WPPRECT_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NECPSST':
            ### Calculate indices for each ensemble member of SPEAR
            NECPSST_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                NECPSST_spear[ens,:] = calc_NECPSST(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = NECPSST_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
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
    np.savez(directorydata + 'ClimateIndices-Tropical_coupledSPEAR-NOAER_%s_1979-2019_trendIncluded.npz' % monq,
                 spear=data_spear,indices=indices,indicesvari=indicesvari,years=years)
