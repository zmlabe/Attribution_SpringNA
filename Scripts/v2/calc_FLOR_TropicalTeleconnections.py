"""
Calculate tropical teleconnections for FLOR

Author    : Zachary M. Labe
Date      : 6 December 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/FLOR/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','May','JFM','FM','FMA','MAM']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    detrend = True
    years = np.arange(1979,2020+1,1)
    
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
        FLOR_lat,FLOR_lon,FLOR = FL.read_FLOR('/work/Zachary.Labe/Data/',
                                              indicesvari,monq,4,'nan',30,'satellite')
        
    
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
        if indices == 'NINO34':
            ### Calculate indices for each ensemble member of FLOR
            NINO34_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                NINO34_FLOR[ens,:] = calc_NINO34(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = NINO34_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'WPPRECT':
            ### Calculate indices for each ensemble member of FLOR
            WPPRECT_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                WPPRECT_FLOR[ens,:] = calc_WPPRECT(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = WPPRECT_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NECPSST':
            ### Calculate indices for each ensemble member of FLOR
            NECPSST_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                NECPSST_FLOR[ens,:] = calc_NECPSST(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = NECPSST_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
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
    np.savez(directorydata + 'ClimateIndices-Tropical_FLOR_%s_1979-2020.npz' % monq,
                 FLOR=data_FLOR,indices=indices,indicesvari=indicesvari,years=years)
