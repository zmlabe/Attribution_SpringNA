"""
Calculate teleconnections for FLOR

Author    : Zachary M. Labe
Date      : 1 December 2022
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
    
    indices = ['NPOcustomZ','ABI','NPO','AL','PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
    indicesvari = ['Z500','Z500','SLP','SLP','Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500']
    
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
    def calc_NPO(data,lat,lon,years):
        """ Calculate NPO for any variables
        """
        ############################################################
        ### Calculated North Pacific Index
        latq_slp1 = np.where((lat >= 55) & (lat <= 72.5))[0]
        lonq_slp1 = np.where((lon >= 180) & (lon <= 220))[0]
        slp1_a = data[:,latq_slp1,:]
        slp1_b = slp1_a[:,:,lonq_slp1]
        lon2_slp1,lat2_slp1 = np.meshgrid(lon[lonq_slp1],lat[latq_slp1])
        slp_1final = UT.calc_weightedAve(slp1_b,lat2_slp1)  
        ############################################################
        ### South SLP box
        latq_slp2 = np.where((lat >= 15) & (lat <= 27.5))[0]
        lonq_slp2 = np.where((lon >= 175) & (lon <= 212.5))[0]
        slp2_a = data[:,latq_slp2,:]
        slp2_b = slp2_a[:,:,lonq_slp2]
        lon2_slp2,lat2_slp2 = np.meshgrid(lon[lonq_slp2],lat[latq_slp2])
        slp_2final = UT.calc_weightedAve(slp2_b,lat2_slp2)
        
        ############################################################
        ### Calculate NPO
        NPO_raw = slp_1final - slp_2final
    
        ############################################################
        ### Normalize NPO by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(NPO_raw[yearq],axis=0)
        mean = np.mean(NPO_raw[yearq],axis=0)
        NPO = (NPO_raw - mean)/std
        
        print('Calculated -----> NPO Index!')
        return NPO
    ###############################################################################
    ###############################################################################
    ###############################################################################
    def calc_NPOcustomZ(data,lat,lon,years):
        """ Calculate NPOcustomZ for any variables
        """
        ############################################################
        ### Calculated North Pacific Index custom Z500
        latq_slp1 = np.where((lat >= 55) & (lat <= 70))[0]
        lonq_slp1 = np.where((lon >= 175) & (lon <= 205))[0]
        slp1_a = data[:,latq_slp1,:]
        slp1_b = slp1_a[:,:,lonq_slp1]
        lon2_slp1,lat2_slp1 = np.meshgrid(lon[lonq_slp1],lat[latq_slp1])
        slp_1final = UT.calc_weightedAve(slp1_b,lat2_slp1)  
        ############################################################
        ### South SLP box
        latq_slp2 = np.where((lat >= 22) & (lat <= 35))[0]
        lonq_slp2 = np.where((lon >= 183) & (lon <= 213))[0]
        slp2_a = data[:,latq_slp2,:]
        slp2_b = slp2_a[:,:,lonq_slp2]
        lon2_slp2,lat2_slp2 = np.meshgrid(lon[lonq_slp2],lat[latq_slp2])
        slp_2final = UT.calc_weightedAve(slp2_b,lat2_slp2)
        
        ############################################################
        ### Calculate NPO
        NPOcustomZ_raw = slp_1final - slp_2final
    
        ############################################################
        ### Normalize NPO by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(NPOcustomZ_raw[yearq],axis=0)
        mean = np.mean(NPOcustomZ_raw[yearq],axis=0)
        NPOcustomZ = (NPOcustomZ_raw - mean)/std
        
        print('Calculated -----> NPOcustomZ Index!')
        return NPOcustomZ
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
    def calc_ABI(data,lat,lon,years):
        """ Calculate Alaska Blocking Index for any variables
        """
        ### Calculate ABI
        lonq = np.where((lon >=180) & (lon <=235))[0]
        latq = np.where((lat >=55) & (lat <=75))[0]
        anomlon = data[:,:,lonq]
        anoms = anomlon[:,latq,:]
        lon2,lat2 = np.meshgrid(lon[lonq],lat[latq])
        ABI_raw = UT.calc_weightedAve(anoms,lat2)
    
        ############################################################
        ### Normalize ABI by 1981-2010
        yearq = np.where((years >= 1981) & (years <= 2010))[0]
        std = np.std(ABI_raw[yearq],axis=0)
        mean = np.mean(ABI_raw[yearq],axis=0)
        ABI = (ABI_raw - mean)/std
        
        print('Calculated -----> ABI Index!')
        return ABI
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
        if indices == 'PNA':
            ### Calculate indices for each ensemble member of FLOR
            pna_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                pna_FLOR[ens,:] = calc_PNA(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = pna_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NAO':
            ### Calculate indices for each ensemble member of FLOR
            nao_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                nao_FLOR[ens,:] = calc_NAO(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = nao_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif any([indices=='PCHT2M',indices=='PCHZ1000',indices=='PCHZ100',indices=='PCHZ50',indices=='PCHZ30']):
            ### Calculate indices for each ensemble member of FLOR
            pch_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                pch_FLOR[ens,:] = calc_PCH(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years)
                
            ### Return common index
            indexz_FLOR = pch_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'UBI':
            ### Calculate indices for each ensemble member of FLOR
            ubi_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                ubi_FLOR[ens,:] = calc_UBI(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
        
            ### Return common index
            indexz_FLOR = ubi_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NPI':
            ### Calculate indices for each ensemble member of FLOR
            npi_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                npi_FLOR[ens,:] = calc_NPI(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
        
            ### Return common index
            indexz_FLOR = npi_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NPO':
            
            ### Calculate indices for each ensemble member of FLOR
            NPO_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                NPO_FLOR[ens,:] = calc_NPO(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
        
            ### Return common index
            indexz_FLOR = NPO_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NPOcustomZ':
            
            ### Calculate indices for each ensemble member of FLOR
            NPOcustomZ_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                NPOcustomZ_FLOR[ens,:] = calc_NPOcustomZ(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
        
            ### Return common index
            indexz_FLOR = NPOcustomZ_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'AL':
            ### Calculate indices for each ensemble member of FLOR
            AL_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                AL_FLOR[ens,:] = calc_AL(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
        
            ### Return common index
            indexz_FLOR = AL_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'SHI':    
            ### Calculate indices for each ensemble member of FLOR
            shi_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                shi_FLOR[ens,:] = calc_SHI(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
                
            ### Return common index
            indexz_FLOR = shi_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'ABI':
            
            ### Calculate indices for each ensemble member of FLOR
            ABI_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                ABI_FLOR[ens,:] = calc_ABI(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
                
            ### Return common index
            indexz_FLOR = ABI_FLOR
        ###############################################################################
        ###############################################################################
        ###############################################################################
        elif indices == 'NINO34':     
            ### Calculate indices for each ensemble member of FLOR
            nino34_FLOR = np.empty((len(FLOR_dt),years.shape[0]))
            for ens in range(len(FLOR_dt)):
                nino34_FLOR[ens,:] = calc_NINO34(FLOR_dt[ens,:,:,:],FLOR_lat,FLOR_lon,years) 
                
            ### Return common index
            indexz_FLOR = nino34_FLOR
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
    np.savez(directorydata + 'ClimateIndices_coupledFLOR_%s_1979-2020.npz' % monq,
                 FLOR=data_FLOR,indices=indices,indicesvari=indicesvari,years=years)
