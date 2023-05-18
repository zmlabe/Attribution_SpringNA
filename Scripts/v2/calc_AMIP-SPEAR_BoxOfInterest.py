"""
Calculate T2M_BoxNA teleconnections for the AMIP experiments

Author    : Zachary M. Labe
Date      : 20 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools
import read_FACTS_AMIPS as AM
import calc_DetrendData as DTT

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/'
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','May','JFM','FM','FMA','MAM','JJA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    scenario = ['spear_obs_rf']
    model_spear = ['SPEAR']
    model = [model_spear]
    detrend = True
    modelunravel = list(itertools.chain(*model))
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
    def calcTeleconnections(scenario,model,indicesvari,indices,monq,years):
        lat = []
        lon = []
        amip = []
        for s in range(len(scenario)):
            late = []
            lone = []
            amipe = []
            for m in range(len(model[s])):
                lat1q,lon1q,lev1q,varq,yearsq = AM.read_FACTS_Experi(scenario[s],model[s][m],
                                                              indicesvari,monq,4,
                                                              'nan','surface')       
                yearsq = np.where((yearsq >= years.min()) & (yearsq <= years.max()))[0]
                varq = varq[:,yearsq,:,:]
            
                late.append(lat1q)
                lone.append(lon1q)
                amipe.append(varq)
            lat.append(late)
            lon.append(lone)
            amip.append(amipe)
            
        ### Organize models  
        spear = np.asarray(amip).squeeze()
        spear_lat = np.asarray(lat).squeeze()
        spear_lon = np.asarray(lon).squeeze()
        
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
        data_spear[iii,:,:] = calcTeleconnections(scenario,
                                                                model,
                                                                indicesvari[iii],
                                                                indices[iii],
                                                                monq,
                                                                years)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Save output 
    directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/'
    np.savez(directoryoutput + 'ClimateIndices-T2M_BoxNA_AMIPSPEAR_%s_1979-2020.npz' % monq,
             spear=data_spear,indices=indices,indicesvari=indicesvari,years=years)
