"""
Calculate EOF for NPO using AMIP-SPEAR

Author    : Zachary M. Labe
Date      : 18 November 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import calc_Utilities as UT
import sys
import read_FACTS_AMIPS as AM
import calc_DetrendData as DT
from eofs.standard import Eof
import scipy.stats as sts
import cmocean

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/NPO/' 
directorydata = '/work/Zachary.Labe/Data/' 
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/NPO/'

### Create months to loop through
sliceperiod = 'April'
yearsrange = np.arange(1979,2020+1,1)

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
    
    return anoms

### Parameters
timelength = 15
eofall = []
pcall = []
varianceall = []
for i,yr in enumerate(range(yearsrange.min(),yearsrange.max()-timelength+1,1)):
    print (i,yr)
    variq = 'SLP'
    sliceshape = 3
    slicenan = 'nan'
    addclimo = True
    yearmin = yr
    yearmax = yr+timelength

    ### Read data
    lat1,lon1,lev1,varn,yearsall = AM.read_FACTS_Experi('spear_obs_rf','SPEAR',
                                                        variq,sliceperiod,4,'nan','surface')  
    
    ### Read only 1979-2020
    yearq = np.where((yearsall >= 1979) & (yearsall <= 2020))[0]
    years = yearsall[yearq]
    var = varn[:,yearq,:,:]
    
    ### Calculate anomalies
    anoms = calc_anomalies(years,var)

    ### Detrend data
    vardt = DT.detrendData(anoms,'surface','yearly')
    yearq = np.where((years >= yearmin) & (years <= yearmax))[0]
    vardtpre = vardt[:,yearq,:,:]

    eofallens = []
    pcallens = []
    varianceallens = []
    for ens in range(len(varn)):
        vardt = vardtpre[ens,:,:,:]
        
        ### Calculate North Pacific Box
        lon2,lat2 = np.meshgrid(lon1,lat1)
        lonq = np.where((lon1 >=120) & (lon1 <=240))[0]
        latq = np.where((lat1 >=20) & (lat1 <=60))[0]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        lon2sq = lon2[latq,:]
        lon2s = lon2sq[:,lonq]
        lat1s = lat1[latq]
        lon1s = lon1[lonq]
        
        anomlon = vardt[:,:,lonq]
        anoms = anomlon[:,latq,:]
        
        ### Calculate EOF (NPO is EOF2 of North Pacific)
        coslat = np.cos(np.deg2rad(lat1s)).clip(0.,1.)
        wgts = np.sqrt(coslat)[...,np.newaxis]
        solver = Eof(anoms,weights=wgts)
        
        ### First 2 eofs
        eofpattern = solver.eofs(neofs=2)
        eof1 = eofpattern[0]
        eof2 = eofpattern[1]
        
        ### Calculate and standardize the PCs
        pc = solver.pcs(npcs=2,pcscaling=0) # not scaled
        pc1 = sts.zscore(pc[:,0])
        pc2 = sts.zscore(pc[:,1])
        
        ### Variance explained by eof pattern
        variance_fraction = solver.varianceFraction(neigs=2)
        variance1 = variance_fraction[0]
        variance2 = variance_fraction[1]
        
        ### Make same sign
        if eof2[2,70] > 0:
            eof2 = eof2 * -1
            pc2 = pc2 * -1
        if np.nanmean(eof1) < 0:
            eof1 = eof1 * -1
            pc1 = pc1 * -1
        
        eofallens.append(eof2)
        pcallens.append(pc2)
        varianceallens.append(variance2)
    eofall.append(eofallens)
    pcall.append(pcallens)
    varianceall.append(varianceallens)
      
minall_lat = []
maxall_lat = []
minall_lon = []
maxall_lon = []
for ens in range(len(varn)):
    ### Calculate location of max node for EOF2
    maxlat = []
    maxlon = []
    for d in range(len(eofall)):
        latwhere = np.unravel_index(eofall[d][ens].argmax(),eofall[d][ens].shape)[0]
        lonwhere = np.unravel_index(eofall[d][ens].argmax(),eofall[d][ens].shape)[1]
        latnode = lat1s[latwhere]
        lonnode = lon1s[lonwhere]
        
        maxlat.append(latnode)
        maxlon.append(lonnode)
        
    ### Calculate location of min node for EOF2
    minlat = []
    minlon = []
    for d in range(len(eofall)):
        latwheremin = np.unravel_index(eofall[d][ens].argmin(),eofall[d][ens].shape)[0]
        lonwheremin = np.unravel_index(eofall[d][ens].argmin(),eofall[d][ens].shape)[1]
        latnodemin = lat1s[latwheremin]
        lonnodemin = lon1s[lonwheremin]
        
        minlat.append(latnodemin)
        minlon.append(lonnodemin)
    
    minall_lat.append(minlat)
    maxall_lat.append(maxlat)
    minall_lon.append(minlon)
    maxall_lon.append(maxlon)

### Save data
np.savez(directoryoutput + 'LatLon_NPO_AMIPSPEAR_%s.npz' % sliceperiod,
         minall_lat=minall_lat,maxall_lat=maxall_lat,
         minall_lon=minall_lon,maxall_lon=maxall_lon)
np.savez(directoryoutput + 'EOFinfo_NPO_AMIPSPEAR_%s.npz' % sliceperiod,
         eofall=eofall,pcall=pcall,varianceall=varianceall,
         timelength=timelength,years=years)

### Save lat/lon grid for plotting
np.savez(directoryoutput + 'LatLonGrid_NPO_AMIPSPEAR_%s.npz' % sliceperiod,
         lat2s=lat2s,lon2s=lon2s)
