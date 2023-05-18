"""
Calculate teleconnections for SPEAR-MED and do not detrend

Author    : Zachary M. Labe
Date      : 15 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_SPEAR_MED as SP

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/coupledSPEAR/'

###############################################################################
###############################################################################
###############################################################################
### Parameters
monthsall = ['January','February','March','April','May','JFM','FM','FMA']
for mn in range(len(monthsall)):
    monq = monthsall[mn]
    detrend = False
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
    def calcTeleconnections(indicesvari,indices,monq,years):
        spear_lat,spear_lon,spear = SP.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',
                                              indicesvari,monq,4,'nan',30,'all')
    
        ### Calculate anomalies 
        spear_dt = calc_anomalies(years,spear) 
        
        ###############################################################################
        ###############################################################################
        ###############################################################################
        if indices == 'PNA':
            ### Calculate indices for each ensemble member of SPEAR
            pna_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                pna_spear[ens,:] = calc_PNA(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = pna_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NAO':
            ### Calculate indices for each ensemble member of SPEAR
            nao_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                nao_spear[ens,:] = calc_NAO(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = nao_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif any([indices=='PCHT2M',indices=='PCHZ1000',indices=='PCHZ100',indices=='PCHZ50',indices=='PCHZ30']):
            ### Calculate indices for each ensemble member of SPEAR
            pch_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                pch_spear[ens,:] = calc_PCH(spear_dt[ens,:,:,:],spear_lat,spear_lon,years)
                
            ### Return common index
            indexz_spear = pch_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'UBI':
            ### Calculate indices for each ensemble member of SPEAR
            ubi_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                ubi_spear[ens,:] = calc_UBI(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_spear = ubi_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NPI':
            ### Calculate indices for each ensemble member of SPEAR
            npi_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                npi_spear[ens,:] = calc_NPI(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_spear = npi_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'AL':
            ### Calculate indices for each ensemble member of SPEAR
            AL_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                AL_spear[ens,:] = calc_AL(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
        
            ### Return common index
            indexz_spear = AL_spear
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'SHI':    
            ### Calculate indices for each ensemble member of SPEAR
            shi_spear = np.empty((len(spear_dt),years.shape[0]))
            for ens in range(len(spear_dt)):
                shi_spear[ens,:] = calc_SHI(spear_dt[ens,:,:,:],spear_lat,spear_lon,years) 
                
            ### Return common index
            indexz_spear = shi_spear
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
    np.savez(directorydata + 'ClimateIndices_coupledSPEAR_%s_1979-2019_trendIncluded.npz' % monq,
                 spear=data_spear,indices=indices,indicesvari=indicesvari,years=years)
