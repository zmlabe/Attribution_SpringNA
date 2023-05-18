"""
Calculate pattern correlations of different variables
 
Author    : Zachary M. Labe
Date      : 20 June 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools
import read_ERA5_monthly1x1 as ER
from scipy.interpolate import griddata as g

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/validationSPEAR/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','spear']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
experi = ['ECHAM5','ESRL-CAM5','SPEAR','ERA5']
model = [model_1880s_rf,model_obs_rf,model_spear]
modelunravel = list(itertools.chain(*model))
variq = 'T2M'
slicemonth = 'FMA'
trendper = 'all'
slicenan = 'nan'
datareader = True
years = np.arange(1979,2019+1,1)

### Read in observations
latobs,lonobs,levobs,var = ER.read_ERA5_monthly1x1(variq,'/work/Zachary.Labe/Data/',slicemonth,
                                years,3,True,slicenan,'surface')

###############################################################################
###############################################################################
###############################################################################    
if datareader == True:
    datatrends = np.load(directorydata + 'MapsOfTrends_AMIPSPEAR_AMIPs_%s_%s_%s.npz' % (slicemonth,
                                                                variq,
                                                                trendper),
                          allow_pickle=True)
    ens_gridold = datatrends['trends']
    lat = datatrends['lat']
    lon = datatrends['lon']
    
def regrid(lat11,lon11,lat21,lon21,var,years):
    """
    Interpolated on selected grid. Reads ERA5 in as 4d with 
    [year,month,lat,lon]
    """
    
    lon1,lat1 = np.meshgrid(lon11,lat11)
    lon2,lat2 = np.meshgrid(lon21,lat21)
    
    varn_re = np.reshape(var,(var.shape[0],(lat1.shape[0]*lon1.shape[1])))   
    varn = np.empty((var.shape[0],lat2.shape[0],lon2.shape[1]))
    
    print('Completed: Start regridding process:')
    for i in range(years.shape[0]):
        z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,:],(lat2,lon2),method='linear')
        varn[i,:,:] = z
        print('Completed: Year %s Regridding---' % (years[i]))
    return varn

### Prepare experiment for calculations
alllats = lat[1]
lat_ech = np.asarray(alllats[0])
lat_cam = np.asarray(alllats[1])
lat_spe = np.asarray(lat[2][0])
alllons = lon[1]
lon_ech = np.asarray(alllons[0])
lon_cam = np.asarray(alllons[1])
lon_spe = np.asarray(lon[2][0])

### Regrid to match climate model data
obs_ech = regrid(latobs,lonobs,lat_ech,lon_ech,var,years)
obs_cam = regrid(latobs,lonobs,lat_cam,lon_cam,var,years)
obs_spe = regrid(latobs,lonobs,lat_spe,lon_spe,var,years)

### Save regridded obs
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/Obs_Regrid/'
np.savez(directoryoutput + 'Regrid_ERA5_%s_%s_ECHAM5.npz' % (variq,slicemonth),data=obs_ech,lat=lat_ech,lon=lon_ech)
np.savez(directoryoutput + 'Regrid_ERA5_%s_%s_ESRL-CAM5.npz' % (variq,slicemonth),data=obs_cam,lat=lat_cam,lon=lon_cam)
np.savez(directoryoutput + 'Regrid_ERA5_%s_%s_SPEAR.npz' % (variq,slicemonth),data=obs_spe,lat=lat_spe,lon=lon_spe)

### Calculate decadal linear trend
yearmin = 1979
yearmax = 2019
trendobs_ech = UT.linearTrendR(obs_ech,years,'surface',yearmin,yearmax)*10
trendobs_cam = UT.linearTrendR(obs_cam,years,'surface',yearmin,yearmax)*10
trendobs_spe = UT.linearTrendR(obs_spe,years,'surface',yearmin,yearmax)*10

### Save regridded trend
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/Obs_Regrid/'
np.savez(directoryoutput + 'Regrid_ERA5_%s-trends_%s_ECHAM5.npz' % (variq,slicemonth),data=trendobs_ech,lat=lat_ech,lon=lon_ech)
np.savez(directoryoutput + 'Regrid_ERA5_%s-trends_%s_ESRL-CAM5.npz' % (variq,slicemonth),data=trendobs_cam,lat=lat_cam,lon=lon_cam)
np.savez(directoryoutput + 'Regrid_ERA5_%s-trends_%s_SPEAR.npz' % (variq,slicemonth),data=trendobs_spe,lat=lat_spe,lon=lon_spe)

### Prepare individual amips to compare
echam = np.asarray(ens_gridold[1][0])
cam = np.asarray(ens_gridold[1][1])
spear = np.asarray(ens_gridold[2][0])

### Function to calculate pattern correlations over different regions
def patternCorrelationRegion(dataobs,datamodel,lat1,lon1,region):
    if region != 'global':
        if region == 'focusarea':
            latmin = 43
            latmax = 60
            lonmin = 240
            lonmax = 295
        elif region == 'arctic':
            latmin = 65
            latmax = 90
            lonmin = 0
            lonmax = 360
        elif region == 'NH':
            latmin = 0
            latmax = 90
            lonmin = 0
            lonmax = 360
        elif region == 'NHextra':
            latmin = 30
            latmax = 90
            lonmin = 0
            lonmax = 360
        elif region == 'NA':
            latmin = 10
            latmax = 80
            lonmin = 180
            lonmax = 300
        elif region == 'US':
            latmin = 24
            latmax = 60
            lonmin = 235
            lonmax = 290
        
        latq = np.where((lat1 > latmin) & (lat1 < latmax))[0]
        latnew = lat1[latq]
        dataobs_new1 = dataobs[latq,:]
        datamodels_new1 = datamodel[latq,:]
    
        if np.max(lon1) < 200:
            print(ValueError('Longitude grid is not correct!'))
            sys.exit()
        lonq = np.where((lon1 > lonmin) & (lon1 < lonmax))[0]
        lonnew = lon1[lonq]
        dataobs_new2 = dataobs_new1[:,lonq]
        datamodels_new2 = datamodels_new1[:,lonq]

        ### Prepare for pattern correlations
        dataobs_new = dataobs_new2
        datamodels_new = datamodels_new2
    
    else:
        latnew = lat1
        lonnew = lon1
        dataobs_new = dataobs
        datamodels_new = datamodel
    
    ### Calculate pattern correlation
    corr = UT.calc_spatialCorr(dataobs_new,datamodels_new,latnew,lonnew,'yesnan')
    
    return corr,latnew,lonnew

### Calculate pattern correlations for different regions
regions = ['global','focusarea','arctic','NH','NHextra','NA','US']
corr_ech_all = []
corr_cam_all = []
corr_spe_all = []
for rr in range(len(regions)):
    
    corr_ech = np.empty((echam.shape[0]))
    for e in range(echam.shape[0]):
        corr_ech[e],latnew,lonnew = patternCorrelationRegion(trendobs_ech,echam[e],lat_ech,lon_ech,regions[rr])
    corr_ech_all.append(corr_ech)
    
    corr_cam = np.empty((cam.shape[0]))
    for e in range(cam.shape[0]):
        corr_cam[e],latnew_cam,lonnew_cam = patternCorrelationRegion(trendobs_cam,cam[e],lat_cam,lon_cam,regions[rr])
    corr_cam_all.append(corr_cam)
    
    corr_spe = np.empty((spear.shape[0]))
    for e in range(spear.shape[0]):
        corr_spe[e],latnew_spe,lonnew_spe = patternCorrelationRegion(trendobs_spe,spear[e],lat_spe,lon_spe,regions[rr])
    corr_spe_all.append(corr_spe)
    
### Save correlations
directorycorr = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/PatternCorrelations/'
np.savez(directorycorr + 'PatternCorrelations_Regional_%s_%s_AMIPSPEAR.npz' % (variq,slicemonth),
         echam=corr_ech_all,cam=corr_cam_all,spear=corr_spe_all)
