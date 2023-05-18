"""
Compare ensemble members with similar trends or even colder
Author    : Zachary M. Labe
Date      : 1 July 2022
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
monthsall = ['January','February','March','April','May','JFM','FM','FMA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    scenario = ['amip_obs_rf','spear']
    model_obs_rf = ['ECHAM5','ESRL-CAM5']
    model_spear = ['spear']
    model = [model_obs_rf,model_spear]
    detrend = True
    modelunravel = list(itertools.chain(*model))
    years = np.arange(1979,2019+1,1)

    indices = ['AL','PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
    indicesvari = ['SLP','Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500']
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Functions for this script
    def calc_PNA(data,lat,lon,years):
        """ Calculate PNA index
        """
        ############################################################
        ### Modified Pointwise Method for the PNA Index from the CPC
        latq_z1 = np.where((lat >= 15) & (lat <=25))[0]
        lonq_z1 = np.where((lon >= 180) & (lon <= 220))[0]
        z1_a = data[:,latq_z1,:]
        z1_b = z1_a[:,:,lonq_z1]
        lon2_z1,lat2_z1 = np.meshgrid(lon[lonq_z1],lat[latq_z1])
        z1_m = UT.calc_weightedAve(z1_b,lat2_z1)  
        ############################################################
        latq_z2 = np.where((lat >= 40) & (lat <=50))[0]
        lonq_z2 = np.where((lon >= 180) & (lon <= 220))[0]
        z2_a = data[:,latq_z2,:]
        z2_b = z2_a[:,:,lonq_z2]
        lon2_z2,lat2_z2 = np.meshgrid(lon[lonq_z2],lat[latq_z2])
        z2_m = UT.calc_weightedAve(z2_b,lat2_z2)
        ############################################################
        latq_z3 = np.where((lat >= 45) & (lat <=60))[0]
        lonq_z3 = np.where((lon >= 235) & (lon <= 255))[0]
        z3_a = data[:,latq_z3,:]
        z3_b = z3_a[:,:,lonq_z3]
        lon2_z3,lat2_z3 = np.meshgrid(lon[lonq_z3],lat[latq_z3])
        z3_m = UT.calc_weightedAve(z3_b,lat2_z3)
        ############################################################
        latq_z4 = np.where((lat >= 25) & (lat <=35))[0]
        lonq_z4 = np.where((lon >= 270) & (lon <= 290))[0]
        z4_a = data[:,latq_z4,:]
        z4_b = z4_a[:,:,lonq_z4]
        lon2_z4,lat2_z4 = np.meshgrid(lon[lonq_z4],lat[latq_z4])
        z4_m = UT.calc_weightedAve(z4_b,lat2_z4)
        
        ############################################################
        ### Calculate PNA
        pna_raw = z1_m - z2_m + z3_m - z4_m
        
        ############################################################
        ### Normalize PNA by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(pna_raw[yearq],axis=0)
        mean = np.mean(pna_raw[yearq],axis=0)
        pna = (pna_raw - mean)/std
        
        print('Calculated -----> PNA Index!')
        return pna 
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_NAO(data,lat,lon,years):
        """ Calculate NAO index
        """
        ############################################################
        ### McKenna et al. 2021, GRL
        lonq1_z1 = np.where((lon >=0) & (lon <=60))[0]
        lonq2_z1 = np.where((lon >= 270) & (lon <= 360))[0]
        lonq_z1 = np.append(lonq1_z1,lonq2_z1)
        latq_z1 = np.where((lat >=20) & (lat <=55))[0]
        z1_a = data[:,latq_z1,:]
        z1_b = z1_a[:,:,lonq_z1]
        lon2_z1,lat2_z1 = np.meshgrid(lon[lonq_z1],lat[latq_z1])
        z1_m = UT.calc_weightedAve(z1_b,lat2_z1)  
        ############################################################
        lonq1_z2 = np.where((lon >=0) & (lon <=60))[0]
        lonq2_z2 = np.where((lon >= 270) & (lon <= 360))[0]
        lonq_z2 = np.append(lonq1_z2,lonq2_z2)
        latq_z2 = np.where((lat >=55) & (lat <=90))[0]
        z2_a = data[:,latq_z2,:]
        z2_b = z2_a[:,:,lonq_z2]
        lon2_z2,lat2_z2 = np.meshgrid(lon[lonq_z2],lat[latq_z2])
        z2_m = UT.calc_weightedAve(z2_b,lat2_z2)
        
        ############################################################
        ### Calculate NAO
        nao_raw = z1_m - z2_m
        
        ############################################################
        ### Normalize NAO by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(nao_raw[yearq],axis=0)
        mean = np.mean(nao_raw[yearq],axis=0)
        nao = (nao_raw - mean)/std
        
        print('Calculated -----> NAO Index!')
        return nao
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_PCH(data,lat,lon,years):
        """ Calculate PCH for any variables
        """
        ############################################################
        ### Calculated polar cap height for >65N
        latq = np.where((lat >= 65) & (lat <=90))[0]
        lonq = np.where((lon >= 0) & (lon <= 360))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        pch_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize PCHB by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(pch_raw[yearq],axis=0)
        mean = np.mean(pch_raw[yearq],axis=0)
        pch = (pch_raw - mean)/std
        
        print('Calculated -----> PCH Index!')
        return pch 
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_NPI(data,lat,lon,years):
        """ Calculate NPI for any variables
        """
        ############################################################
        ### Calculated North Pacific Index
        lonq = np.where((lon >=160) & (lon <=220))[0]
        latq = np.where((lat >=30) & (lat <=65))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        npi_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize NPI by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(npi_raw[yearq],axis=0)
        mean = np.mean(npi_raw[yearq],axis=0)
        npi = (npi_raw - mean)/std
        
        print('Calculated -----> NPI Index!')
        return npi 
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_AL(data,lat,lon,years):
        """ Calculate AL for any variables
        """
        ############################################################
        ### Calculated Aleutian Low Index
        lonq = np.where((lon >=160) & (lon <=200))[0]
        latq = np.where((lat >=45) & (lat <=65))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        AL_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize AL by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(AL_raw[yearq],axis=0)
        mean = np.mean(AL_raw[yearq],axis=0)
        AL = (AL_raw - mean)/std
        
        print('Calculated -----> AL Index!')
        return AL
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_NINO34(data,lat,lon,years):
        """ Calculate NINO34 for any variables
        """
        ############################################################
        ### Calculated NINO3.4 Index
        lonq = np.where((lon >=190) & (lon <=240))[0]
        latq = np.where((lat >=-5) & (lat <=5))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        nino34_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize NINO34 by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(nino34_raw[yearq],axis=0)
        mean = np.mean(nino34_raw[yearq],axis=0)
        nino34 = (nino34_raw - mean)/std
        
        print('Calculated -----> NINO34 Index!')
        return nino34 
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_SHI(data,lat,lon,years):
        """ Calculate SHI for any variables
        """
        ############################################################
        ### Calculated Siberian High Index
        lonq = np.where((lon >=80) & (lon <=120))[0]
        latq = np.where((lat >=40) & (lat <=65))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        shi_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize Siberian High Index by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(shi_raw[yearq],axis=0)
        mean = np.mean(shi_raw[yearq],axis=0)
        shi = (shi_raw - mean)/std
        
        print('Calculated -----> SHI Index!')
        return shi 
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_UBI(data,lat,lon,years):
        """ Calculate Ural Blocking Index for any variables
        """
        ############################################################
        ### Calculated Ural Blocking Index
        lonq1 = np.where((lon >=0) & (lon <=90))[0]
        lonq2 = np.where((lon >= 330) & (lon <= 360))[0]
        lonq = np.append(lonq1,lonq2)
        latq = np.where((lat >=45) & (lat <=80))[0]
        var_a = data[:,latq,:]
        var_b = var_a[:,:,lonq]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        ubi_raw = UT.calc_weightedAve(var_b,lat2)  
    
        ############################################################
        ### Normalize UBI by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(ubi_raw[yearq],axis=0)
        mean = np.mean(ubi_raw[yearq],axis=0)
        ubi = (ubi_raw - mean)/std
        
        print('Calculated -----> UBI Index!')
        return ubi
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
        echam = amip[0][0]
        echam_lat = lat[0][0]
        echam_lon = lon[0][0]
        
        cam = amip[0][1]
        cam_lat = lat[0][1]
        cam_lon = lon[0][1]
        
        spear = amip[1][0]
        spear_lat = lat[1][0]
        spear_lon = lon[1][0]
        
        ### Calculate anomalies 
        echam_anom = calc_anomalies(years,echam) 
        cam_anom = calc_anomalies(years,cam) 
        spear_anom = calc_anomalies(years,spear) 
        
        ### Detrend data
        if detrend == True:
            echam_dt = DTT.detrendData(echam_anom,'surface','yearly')
            cam_dt = DTT.detrendData(cam_anom,'surface','yearly')
            spear_dt = DTT.detrendData(spear_anom,'surface','yearly')
        else:
            echam_dt = echam_anom
            cam_dt = cam_anom
            spear_dt = spear_anom
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'PNA':
            ### Calculate indices for each ensemble member of ECHAM5
            pna_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                pna_echam[ens,:] = calc_PNA(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            pna_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                pna_cam[ens,:] = calc_PNA(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            pna_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                pna_spear[ens,:] = calc_PNA(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_echam = pna_echam
            indexz_cam = pna_cam
            indexz_spear = pna_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NAO':
            ### Calculate indices for each ensemble member of ECHAM5
            nao_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                nao_echam[ens,:] = calc_NAO(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            nao_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                nao_cam[ens,:] = calc_NAO(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            nao_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                nao_spear[ens,:] = calc_NAO(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_echam = nao_echam
            indexz_cam = nao_cam
            indexz_spear = nao_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif any([indices=='PCHT2M',indices=='PCHZ1000',indices=='PCHZ100',indices=='PCHZ50',indices=='PCHZ30']):
            ### Calculate indices for each ensemble member of ECHAM5
            pch_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                pch_echam[ens,:] = calc_PCH(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            pch_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                pch_cam[ens,:] = calc_PCH(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            pch_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                pch_spear[ens,:] = calc_PCH(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_echam = pch_echam
            indexz_cam = pch_cam
            indexz_spear = pch_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'UBI':
            ### Calculate indices for each ensemble member of ECHAM5
            ubi_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                ubi_echam[ens,:] = calc_UBI(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            ubi_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                ubi_cam[ens,:] = calc_UBI(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            ubi_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                ubi_spear[ens,:] = calc_UBI(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_echam = ubi_echam
            indexz_cam = ubi_cam
            indexz_spear = ubi_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'AL':
            ### Calculate indices for each ensemble member of ECHAM5
            AL_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                AL_echam[ens,:] = calc_AL(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            AL_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                AL_cam[ens,:] = calc_AL(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            AL_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                AL_spear[ens,:] = calc_AL(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_echam = AL_echam
            indexz_cam = AL_cam
            indexz_spear = AL_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NPI':
            ### Calculate indices for each ensemble member of ECHAM5
            npi_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                npi_echam[ens,:] = calc_NPI(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            npi_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                npi_cam[ens,:] = calc_NPI(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            npi_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                npi_spear[ens,:] = calc_NPI(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_echam = npi_echam
            indexz_cam = npi_cam
            indexz_spear = npi_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'SHI':
            ### Calculate indices for each ensemble member of ECHAM5
            shi_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                shi_echam[ens,:] = calc_SHI(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            shi_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                shi_cam[ens,:] = calc_SHI(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            shi_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                shi_spear[ens,:] = calc_SHI(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
                
            ### Return common index
            indexz_echam = shi_echam
            indexz_cam = shi_cam
            indexz_spear = shi_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NINO34':
            ### Calculate indices for each ensemble member of ECHAM5
            nino34_echam = np.empty((len(echam_dt),years.shape[0]))
            for ens in range(len(echam_dt)):
                nino34_echam[ens,:] = calc_NINO34(echam_dt[ens,:,:,:],echam_lat,echam_lon,years)
                
            ### Calculate indices for each ensemble member of CAM5
            nino34_cam = np.empty((len(cam_dt),years.shape[0]))
            for ens in range(len(cam_dt)):
                nino34_cam[ens,:] = calc_NINO34(cam_dt[ens,:,:,:],cam_lat,cam_lon,years)
            
            ### Calculate indices for each ensemble member of SPEAR
            nino34_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                nino34_spear[ens,:] = calc_NINO34(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
                
            ### Return common index
            indexz_echam = nino34_echam
            indexz_cam = nino34_cam
            indexz_spear = nino34_spear
        else:
            print(ValueError('SOMETHING IS WRONG WITH THE INDEX SELECTED!'))
            sys.exit()
        return indexz_echam,indexz_cam,indexz_spear
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Begin loop for indices for calculations
    if scenario[0] == 'amip_obs_rf':
        data_echam = np.empty((len(indices),50,years.shape[0]))
        data_cam = np.empty((len(indices),40,years.shape[0]))
        data_spear = np.empty((len(indices),3,years.shape[0]))
    elif scenario[0] == 'amip_1880s_rf':
        data_echam = np.empty((len(indices),50,years.shape[0]))
        data_cam = np.empty((len(indices),40,years.shape[0]))
        data_spear = np.empty((len(indices),3,years.shape[0]))
    for iii in range(len(indices)):
        data_echam[iii,:,:],data_cam[iii,:,:],data_spear[iii,:,:] = calcTeleconnections(scenario,
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
        np.savez(directoryoutput + 'ClimateIndices_ValidationSPEAR-1880rf_%s_1979-2019.npz' % monq,
                 echam=data_echam,cam=data_cam,spear=data_spear,
                 indices=indices,indicesvari=indicesvari,years=years)
    else:
        np.savez(directoryoutput + 'ClimateIndices_ValidationSPEAR_%s_1979-2019.npz' % monq,
                 echam=data_echam,cam=data_cam,spear=data_spear,
                 indices=indices,indicesvari=indicesvari,years=years)
