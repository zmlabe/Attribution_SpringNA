"""
Plot histogram of trends over different time periods compared to SMILEs CMIP6

Author    : Zachary M. Labe
Date      : 8 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_MIROC6_LE_LOWS as LE
import read_SMHI_LE_LOWS as LE2
import read_MPI_ESM12_LE as FL
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/MMLEA/'

### Parameters
variq = 'T2M'
sliceperiod = 'April'
yearssat = np.arange(1979,2020+1,1)
slicenan = 'nan'
addclimo = True
newon = True
datareader = True

### Read data
if datareader == True:
    lat1l,lon1l,lens1 = LE.read_MIROC6_LE_LOWS('/work/Zachary.Labe/Data/',variq,sliceperiod,
                                      4,slicenan,50,'satellite')
    lat1l2,lon1l2,lens2 = LE2.read_SMHI_LE_LOWS('/work/Zachary.Labe/Data/',variq,sliceperiod,
                                       4,slicenan,50,'satellite')
    lat1s,lon1s,flor = FL.read_MPI_ESM12_LE('/work/Zachary.Labe/Data/',variq,sliceperiod,
                                        4,slicenan,30,'satellite')
  
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

def calc_regionalAve(regionn,lat1,lon1,data):
    if regionn == 'CANMID':
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
        
    if data.ndim == 3:
        meanlat = data[:,lat1q,:]
        meanbox = meanlat[:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    elif data.ndim == 4:
        meanlat = data[:,:,lat1q,:]
        meanbox = meanlat[:,:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
        
    ### Calculate timeseries
    mean = UT.calc_weightedAve(meanbox,lat2q)

    return mean,lat1a,lat2q,lon1a,lon2q

def calcTrends(years,yrmin,yrmax,data):
    yearq = np.where((years >= yrmin) & (years <= yrmax))[0]
    yearTrend = years[yearq]
    print(yearTrend)
    
    if data.ndim == 1:
        dataTrend = data[yearq]
        slope,intercept,r,p,se = sts.linregress(yearTrend,dataTrend)
        trendline = slope*yearTrend + intercept
    elif data.ndim == 2:
        dataTrend = data[:,yearq]
        
        slope = np.empty((data.shape[0]))
        intercept = np.empty((data.shape[0]))
        r = np.empty((data.shape[0]))
        p = np.empty((data.shape[0]))
        se = np.empty((data.shape[0]))
        trendline = np.empty((data.shape[0],len(yearTrend)))
        for i in range(data.shape[0]):
            slopeq,interceptq,rq,pq,seq = sts.linregress(yearTrend,dataTrend[i,:])
            slope[i] = slopeq
            intercept[i] = interceptq
            r[i] = rq
            p[i] = pq
            se[i] = seq
            trendline[i,:] = slopeq*yearTrend + interceptq
    else:
        print(ValueError('WRONG DIMENSIONS OF TREND INPUT!!!'))
        
    return slope,intercept,r,p,se,trendline,yearTrend

lens1_anom = calc_anomalies(yearssat,lens1)
lens2_anom = calc_anomalies(yearssat,lens2)
flor_anom = calc_anomalies(yearssat,flor)

lens_reg1,latr1l,latr2l,lonr1l,lonr2l = calc_regionalAve('CANMID',lat1l,lon1l,lens1_anom)
lens_reg2,latr1l2,latr2l2,lonr1l2,lonr2l2 = calc_regionalAve('CANMID',lat1l2,lon1l2,lens2_anom)
flor_reg,latrs1,latrs2,lonrs1,lonrs2 = calc_regionalAve('CANMID',lat1s,lon1s,flor_anom)

yrmin3 = 1979
yrmax3 = 2020
slope_lens3,intercept_lens3,r_lens3,p_lens3,se_lens3,trendline_lens3,yearTrend_lens3 = calcTrends(yearssat,yrmin3,yrmax3,lens_reg1)
slope_lens32,intercept_lens32,r_lens32,p_lens32,se_lens32,trendline_lens32,yearTrend_lens32 = calcTrends(yearssat,yrmin3,yrmax3,lens_reg2)
slope_flor3,intercept_flor3,r_flor3,p_flor3,se_flor3,trendline_flor3,yearTrend_flor3 = calcTrends(yearssat,yrmin3,yrmax3,flor_reg)

directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/MMLEA/'
np.savetxt(directoryoutput + 'Slopes_%s_MIRO6_LE_LOWS_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_lens3)
np.savetxt(directoryoutput + 'Slopes_%s_SMHI_LE_LOWS_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_lens32)
np.savetxt(directoryoutput + 'Slopes_%s_MPI_ESM12_LE_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_flor3)

