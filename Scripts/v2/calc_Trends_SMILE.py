"""
Plot histogram of trends over different time periods compared to SMILEs

Author    : Zachary M. Labe
Date      : 23 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
import read_SMILE_historical as SM
import read_LENS1 as LE
import read_LENS2 as LE2
import read_SPEAR_MED as SP
import read_SPEAR_MED_NOAER as AER
import read_SPEAR_MED_NATURAL as NAT
import read_FLOR as FL
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'FMA'
years = np.arange(1979,2021+1,1)
yearssat = np.arange(1979,2020+1,1)
slicenan = 'nan'
addclimo = True
newon = True
datareader = True

### Read data
if datareader == True:
    lat1e,lon1e,lev1,obs = ERA.read_ERA5_monthly1x1(variq,'/work/Zachary.Labe/Data/',sliceperiod,
                                             years,3,addclimo,slicenan,'surface')
    lat1,lon1,gfdlc = SM.read_SMILEhistorical('/work/Zachary.Labe/Data/SMILE/','GFDL_CM3',variq,
                                                sliceperiod,4, slicenan,20)
    lat1,lon1,gfdlem = SM.read_SMILEhistorical('/work/Zachary.Labe/Data/SMILE/','GFDL_ESM2M',variq,
                                                sliceperiod,4,slicenan,20)
    lat1l,lon1l,lens1 = LE.read_LENS1('/work/Zachary.Labe/Data/',variq,sliceperiod,
                                      4,slicenan,40,'satellite')
    lat1l2,lon1l2,lens2 = LE2.read_LENS2('/work/Zachary.Labe/Data/',variq,sliceperiod,
                                       4,slicenan,100,'satellite')
    lat1s,lon1s,flor = FL.read_FLOR('/work/Zachary.Labe/Data/',variq,sliceperiod,
                                        4,slicenan,30,'satellite')
    lat1s,lon1s,spear = SP.read_SPEAR_MED('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED/monthly/',variq,
                                            sliceperiod,4,slicenan,30,'all')
    lat1sa,lon1sa,spearnoaer = AER.read_SPEAR_MED_NOAER('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED_NOAER/monthly/',variq,
                                            sliceperiod,4,slicenan,12,'satellite')
    lat1snat,lon1snat,spearnat = NAT.read_SPEAR_MED_NATURAL('/work/Zachary.Labe/Data/SPEAR/SPEAR_MED_NATURAL/monthly/',variq,
                                            sliceperiod,4,slicenan,30,'satellite')
  
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

obs_anom = calc_anomalies(years,obs)
gfdlc_anom = calc_anomalies(years,gfdlc)
gfdlem_anom = calc_anomalies(years,gfdlem)
lens1_anom = calc_anomalies(yearssat,lens1)
lens2_anom = calc_anomalies(yearssat,lens2)
spear_anom = calc_anomalies(yearssat,spear)
spearnoaer_anom = calc_anomalies(yearssat,spearnoaer)
spearnat_anom = calc_anomalies(yearssat,spearnat)
flor_anom = calc_anomalies(yearssat,flor)

obs_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1e,lon1e,obs_anom)
gfdlc_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,gfdlc_anom)
gfdlem_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,gfdlem_anom)
lens_reg1,latr1l,latr2l,lonr1l,lonr2l = calc_regionalAve('CANMID',lat1l,lon1l,lens1_anom)
lens_reg2,latr1l2,latr2l2,lonr1l2,lonr2l2 = calc_regionalAve('CANMID',lat1l2,lon1l2,lens2_anom)
spear_reg,latrs1,latrs2,lonrs1,lonrs2 = calc_regionalAve('CANMID',lat1s,lon1s,spear_anom)
spearnoaer_reg,latrsa1,latrsa2,lonrsa1,lonrsa2 = calc_regionalAve('CANMID',lat1sa,lon1sa,spearnoaer_anom)
spearnat_reg,latrsnat1,latrsnat2,lonrsnat1,lonrsnat2 = calc_regionalAve('CANMID',lat1snat,lon1snat,spearnat_anom)
flor_reg,latrs1,latrs2,lonrs1,lonrs2 = calc_regionalAve('CANMID',lat1s,lon1s,flor_anom)

yrmin3 = 1979
yrmax3 = 2020
slope_obs3,intercept_obs3,r_obs3,p_obs3,se_obs3,trendline_obs3,yearTrend_obs3 = calcTrends(years,yrmin3,yrmax3,obs_reg)
slope_gfdlc3,intercept_gfdlc3,r_gfdlc3,p_gfdlc3,se_gfdlc3,trendline_gfdlc3,yearTrend_gfdlc3 = calcTrends(yearssat,yrmin3,yrmax3,gfdlc_reg)
slope_gfdlem3,intercept_gfdlem3,r_gfdlem3,p_gfdlem3,se_gfdlem3,trendline_gfdlem3,yearTrend_gfdlem3 = calcTrends(yearssat,yrmin3,yrmax3,gfdlem_reg)
slope_lens3,intercept_lens3,r_lens3,p_lens3,se_lens3,trendline_lens3,yearTrend_lens3 = calcTrends(yearssat,yrmin3,yrmax3,lens_reg1)
slope_lens32,intercept_lens32,r_lens32,p_lens32,se_lens32,trendline_lens32,yearTrend_lens32 = calcTrends(yearssat,yrmin3,yrmax3,lens_reg2)
slope_spear3,intercept_spear3,r_spear3,p_spear3,se_spear3,trendline_spear3,yearTrend_spear3 = calcTrends(yearssat,yrmin3,yrmax3,spear_reg)
slope_spearnoaer3,intercept_spearnoaer3,r_spearnoaer3,p_spearnoaer3,se_spearnoaer3,trendline_spearnoaer3,yearTrend_spearnoaer3 = calcTrends(yearssat,yrmin3,yrmax3,spearnoaer_reg)
slope_spearnat3,intercept_spearnat3,r_spearnat3,p_spearnat3,se_spearnoaer3,trendline_spearnat3,yearTrend_spearnat3 = calcTrends(yearssat,yrmin3,yrmax3,spearnat_reg)
slope_flor3,intercept_flor3,r_flor3,p_flor3,se_flor3,trendline_flor3,yearTrend_flor3 = calcTrends(yearssat,yrmin3,yrmax3,flor_reg)

directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/MMLEA/'
np.savetxt(directoryoutput + 'Slopes_%s_obs_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),[slope_obs3])
np.savetxt(directoryoutput + 'Slopes_%s_gfdlc_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_gfdlc3)
np.savetxt(directoryoutput + 'Slopes_%s_gfdlem_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_gfdlem3)
np.savetxt(directoryoutput + 'Slopes_%s_lens1_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_lens3)
np.savetxt(directoryoutput + 'Slopes_%s_lens2_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_lens32)
np.savetxt(directoryoutput + 'Slopes_%s_spear_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_spear3)
np.savetxt(directoryoutput + 'Slopes_%s_noaerspear_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_spearnoaer3)
np.savetxt(directoryoutput + 'Slopes_%s_natural_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_spearnat3)
np.savetxt(directoryoutput + 'Slopes_%s_FLOR_%s-%s.txt' % (sliceperiod,yrmin3,yrmax3),slope_flor3)

