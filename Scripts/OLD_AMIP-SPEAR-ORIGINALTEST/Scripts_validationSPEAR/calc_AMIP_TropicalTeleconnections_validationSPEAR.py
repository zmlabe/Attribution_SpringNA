"""
Calculate tropical teleconnections for the SPEAR-AMIP experiments

Author    : Zachary M. Labe
Date      : 14 July 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/validationSPEAR/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','JFM','FM','FMA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    scenario = ['spear']
    model_spear = ['spear']
    model = [model_spear]
    detrend = True
    modelunravel = list(itertools.chain(*model))
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
        spear = amip[0][0]
        spear_lat = lat[0][0]
        spear_lon = lon[0][0]
        
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
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE INDEX SELECTED!'))
            sys.exit()
        return indexz_spear
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Begin loop for indices for calculations
    data_spear = np.empty((len(indices),3,years.shape[0]))
    for iii in range(len(indices)):
        data_spear[iii,:,:] = calcTeleconnections(scenario,model,indicesvari[iii],
                                                  indices[iii],monq,years)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Save output 
    directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/Teleconnections/'
    np.savez(directoryoutput + 'ClimateIndices-Tropical_SPEARONLY_%s_1979-2019.npz' % monq,
              spear=data_spear,indices=indices,indicesvari=indicesvari,years=years)
