"""
Plot histogram of trends over different time periods compared to SMILEs

Author    : Zachary M. Labe
Date      : 23 May 2022
"""

import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_FACTS_AMIPS as AM
import scipy.stats as sts
import itertools

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','spear']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
model = [model_1880s_rf,model_obs_rf,model_spear]
variq = 'T2M'
slicemonth= 'FMA'
slicenan = 'nan'
addclimo = True
datareader = True
yearAllData = np.arange(1979,2019+1,1)

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
    
    if data.ndim == 1:
        dataTrend = data[yearq]
        slope,intercept,r,p,se = sts.linregress(yearTrend,dataTrend)
    elif data.ndim == 2:
        dataTrend = data[:,yearq]
        slope = np.empty((data.shape[0]))
        intercept = np.empty((data.shape[0]))
        r = np.empty((data.shape[0]))
        p = np.empty((data.shape[0]))
        se = np.empty((data.shape[0]))
        for i in range(data.shape[0]):
            slopeq,interceptq,rq,pq,seq = sts.linregress(yearTrend,dataTrend[i,:])
            slope[i] = slopeq
            intercept[i] = interceptq
            r[i] = rq
            p[i] = pq
            se[i] = seq
    else:
        print(ValueError('WRONG DIMENSIONS OF TREND INPUT!!!'))
        
    return slope,intercept,r,p,se

def trendPeriods(years,trendper,data):
    slope = []
    intercept = []
    r = []
    p = []
    se = []
    for s in range(len(scenario)):
        slopee = []
        intercepte = []
        re = []
        pe = []
        see = []
        for m in range(len(model[s])):
            yearModel = years
            if trendper == 'all':
                yrmin = 1979
                yrmax = yearModel.max()
            elif trendper == 'recent':
                yrmin = 2000
                yrmax = yearModel.max() 
            print(yrmin,yrmax)
            slopeq,interceptq,rq,pq,seq = calcTrends(yearModel,yrmin,yrmax,data[s][m])
            slopee.append(slopeq)
            intercepte.append(interceptq)
            re.append(rq)
            pe.append(pq)
            see.append(seq)
        slope.append(slopee)
        intercept.append(intercepte)
        r.append(re)
        p.append(pe)
        se.append(see)
    return slope,intercept,r,p,se

### Read data
if datareader == True:
    lat = []
    lon = []
    lat2 = []
    lon2 = []
    lat1all = []
    lon1all = []
    meananom = []
    year = []
    ens = []
    lenyear = []
    amip = []
    for s in range(len(scenario)):
        late = []
        lone = []
        lat2e = []
        lon2e = []
        lat1alle = []
        lon1alle = []
        meananome = []
        yeare = []
        ense = []
        lenyeare = []
        amipe = []
        for m in range(len(model[s])):
            lat1q,lon1q,lev1q,varq,yearsq = AM.read_FACTS_Experi(scenario[s],model[s][m],
                                                         variq,slicemonth,4,
                                                         slicenan,'surface')
            
            yearAllDataq = np.where((yearsq >= yearAllData.min()) & (yearsq <= yearAllData.max()))[0]
            varq = varq[:,yearAllDataq,:,:]
            print(yearAllDataq)
            
            lat1alle.append(lat1q)
            lon1alle.append(lon1q)
            
            ### Calculate anomalies
            anomsq = calc_anomalies(yearAllData,varq)
            
            ### Calculate regional mean
            mean,lat1r,lat2r,lon1r,lon2r = calc_regionalAve('CANMID',lat1q,
                                                            lon1q,anomsq)
            
            late.append(lat1r)
            lone.append(lon1r)
            lat2e.append(lat2r)
            lon2e.append(lon2r)
            meananome.append(mean)
            yeare.append(yearsq)
            ense.append(mean.shape[0])
            lenyeare.append(mean.shape[1])
            amipe.append(varq)
        lat.append(late)
        lon.append(lone)
        lat2.append(lat2e)
        lon2.append(lon2e)
        lat1all.append(lat1alle)
        lon1all.append(lon1alle)
        meananom.append(meananome)
        year.append(yeare)
        ens.append(ense)
        lenyear.append(lenyeare)
        amip.append(amipe)

### Calculate trends
slope_all,intercept_all,r_all,p_all,se_all = trendPeriods(yearAllData,'all',meananom)
slope_recent,intercept_recent,r_recent,p_recent,se_recent = trendPeriods(yearAllData,'recent',meananom)

trendsall = [slope_all,slope_recent]
interceptall = [intercept_all,intercept_recent]
rall = [r_all,r_recent]
pall = [p_all,p_recent]
seall = [se_all,se_recent]

# ### Save metadata and trends for AMIPS
# directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/'
# np.savez(directoryoutput + 'AMIP_SPEARval_AllData_%s_%s.npz' % (slicemonth,variq),amip=amip,
#         lat=lat1all,lon=lon1all)
# np.savez(directoryoutput + 'AMIP_SPEARval_Trends_%s_%s.npz' % (slicemonth,variq),trend=trendsall,
#          intercept=interceptall,r=rall,p=pall,se=seall)
