"""
Calculate SND_BoxNA teleconnections for the AMIP experiments

Author    : Zachary M. Labe
Date      : 21 July 2022
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
monthsall = ['January','February','March','April','May','JFM','FM','FMA','MAM']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    scenario = ['amip_obs_rf','spear_obs_rf']
    model_obs_rf = ['ESRL-CAM5']
    model_spear = ['SPEAR']
    model = [model_obs_rf,model_spear]
    detrend = True
    modelunravel = list(itertools.chain(*model))
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
        cam = amip[0][0]
        cam_lat = lat[0][0]
        cam_lon = lon[0][0]
        
        spear = amip[1][0]
        spear_lat = lat[1][0]
        spear_lon = lon[1][0]
        
        ### Calculate anomalies 
        cam_anom = calc_anomalies(years,cam) 
        spear_anom = calc_anomalies(years,spear) 
        
        ### Detrend data
        if detrend == True:
            cam_dt = DTT.detrendData(cam_anom,'surface','yearly')
            spear_dt = DTT.detrendData(spear_anom,'surface','yearly')
        else:
            cam_dt = cam_anom
            spear_dt = spear_anom
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'SND_BoxNA':
            ### Calculate indices for each ensemble member of CAM5
            SND_BoxNA_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                SND_BoxNA_cam[ens,:] = calc_SND_BoxNA(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            SND_BoxNA_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                SND_BoxNA_spear[ens,:] = calc_SND_BoxNA(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_cam = SND_BoxNA_cam
            indexz_spear = SND_BoxNA_spear
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE INDEX SELECTED!'))
            sys.exit()
        return indexz_cam,indexz_spear
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Begin loop for indices for calculations
    if scenario[0] == 'amip_obs_rf':
        data_cam = np.empty((len(indices),40,years.shape[0]))
        data_spear = np.empty((len(indices),30,years.shape[0]))
    elif scenario[0] == 'amip_1880s_rf':
        data_cam = np.empty((len(indices),30,years.shape[0]))
        data_spear = np.empty((len(indices),30,years.shape[0]))
    for iii in range(len(indices)):
        data_cam[iii,:,:],data_spear[iii,:,:] = calcTeleconnections(scenario,
                                                                model,
                                                                indicesvari[iii],
                                                                indices[iii],
                                                                monq,
                                                                years)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Save output 
    directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/Teleconnections/'
    if scenario[0] == 'amip_1880s_rf':
        np.savez(directoryoutput + 'ClimateIndices-SND_BoxNA_AMIPSPEAR-1880rf_%s_1979-2019.npz' % monq,
                 cam=data_cam,spear=data_spear,
                 indices=indices,indicesvari=indicesvari,years=years)
    else:
        np.savez(directoryoutput + 'ClimateIndices-SND_BoxNA_AMIPSPEAR_%s_1979-2019.npz' % monq,
                 cam=data_cam,spear=data_spear,
                 indices=indices,indicesvari=indicesvari,years=years)
