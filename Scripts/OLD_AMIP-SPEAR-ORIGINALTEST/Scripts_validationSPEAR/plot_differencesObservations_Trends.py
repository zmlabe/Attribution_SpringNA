"""
Plot differences in trends from observational datasets
Author    : Zachary M. Labe
Date      : 23 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthly as ERA
import read_NCEP2 as NP
import read_BEST as B
import read_GISTEMP as G
import read_NOAAglobaltemp as NO
import scipy.stats as sts
import sys
import palettable.cubehelix as cm

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'April'
years = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True

### Read data
lat1,lon1,era5 = ERA.read_ERA5_monthly(variq,'/work/Zachary.Labe/Data/ERA5_19x25/',sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,False)
era5 = era5[:,:,:]
lat1,lon1,ncep = NP.read_NCEP2('/work/Zachary.Labe/Data/NCEP2/',sliceperiod,years,sliceshape,addclimo,slicenan,False)
lat1,lon1,giss = G.read_GISTEMP('/work/Zachary.Labe/Data/GISTEMP/',sliceperiod,years,sliceshape,addclimo,slicenan,False)
lat1b,lon1b,best = B.read_BEST('/work/Zachary.Labe/Data/BEST/',sliceperiod,years,sliceshape,addclimo,slicenan)
lat1,lon1,noaa = NO.read_NOAA('/work/Zachary.Labe/Data/NOAAGlobalTemp/',sliceperiod,years,sliceshape,addclimo,slicenan,False)

### Calculate trends
def calcTrend(years,yrmin,yrmax,data,lat1,lon1):
    yearq = np.where((years >= yrmin) & (years <= yrmax))[0]
    yearTrend = years[yearq]
    print(yearTrend)
    
    ### Calculate region
    la1 = 43
    la2 = 60
    lo1 = 240
    lo2 = 295
    lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
    lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
    
    meanlat = data[:,lat1q,:]
    meanbox = meanlat[:,:,lon1q]
    lon1a = lon1[lon1q]
    lat1a = lat1[lat1q]
    lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    
    ### Calculate timeseries
    mean = UT.calc_weightedAve(meanbox,lat2q)
    
    ### Calculate statistics
    slope, intercept, r, p, se = sts.linregress(yearTrend,mean[yearq])
    trendline = slope*years + intercept

    ### Decadal trend
    dectrend = slope*10
    
    return dectrend,intercept,r,p,se,trendline

dectrend_e1,intercept_e1,r_e1,p_e1,se_e1,trendline_e1 = calcTrend(years,1979,2019,era5,lat1,lon1)
dectrend_n1,intercept_n1,r_n1,p_n1,se_n1,trendline_n1 = calcTrend(years,1979,2019,ncep,lat1,lon1)
dectrend_g1,intercept_g1,r_g1,p_g1,se_g1,trendline_g1 = calcTrend(years,1979,2019,giss,lat1,lon1)
dectrend_b1,intercept_b1,r_b1,p_b1,se_b1,trendline_b1 = calcTrend(years,1979,2019,best,lat1b,lon1b)
dectrend_no1,intercept_no1,r_no1,p_no1,se_no1,trendline_no1 = calcTrend(years,1979,2019,noaa,lat1,lon1)

dectrend_e2,intercept_e2,r_e2,p_e2,se_e2,trendline_e2 = calcTrend(years,2000,2019,era5,lat1,lon1)
dectrend_n2,intercept_n2,r_n2,p_n2,se_n2,trendline_n2 = calcTrend(years,2000,2019,ncep,lat1,lon1)
dectrend_g2,intercept_g2,r_g2,p_g2,se_g2,trendline_g2 = calcTrend(years,2000,2019,giss,lat1,lon1)
dectrend_b2,intercept_b2,r_b2,p_b2,se_b2,trendline_b2 = calcTrend(years,2000,2019,best,lat1b,lon1b)
dectrend_no2,intercept_no2,r_no2,p_no2,se_no2,trendline_no2 = calcTrend(years,2000,2019,noaa,lat1,lon1)

###############################################################################
###############################################################################
###############################################################################               
### Plot Figure
### Adjust axes in time series plots 
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',labelsize=6,pad=4)
ax.tick_params(axis='y',labelsize=6,pad=1.5)
ax.yaxis.grid(zorder=1,color='darkgrey',alpha=0.35,clip_on=False)

plt.axhline(0,color='dimgrey',linestyle='-',linewidth=2,zorder=1)

cmap = cm.cubehelix3_16_r.mpl_colormap

plt.scatter(1,dectrend_e1,color=cmap(0.1),s=70,alpha=0.5,edgecolor='dimgrey',zorder=2)
plt.scatter(1,dectrend_e2,color=cmap(0.1),s=70,alpha=1,edgecolor='dimgrey',zorder=2)

plt.scatter(2,dectrend_n1,color=cmap(0.3),s=70,alpha=0.5,edgecolor='dimgrey',zorder=2)
plt.scatter(2,dectrend_n2,color=cmap(0.3),s=70,alpha=1,edgecolor='dimgrey',zorder=2)

plt.scatter(3,dectrend_g1,color=cmap(0.5),s=70,alpha=0.5,edgecolor='dimgrey',zorder=2)
plt.scatter(3,dectrend_g2,color=cmap(0.5),s=70,alpha=1,edgecolor='dimgrey',zorder=2)

plt.scatter(4,dectrend_b1,color=cmap(0.7),s=70,alpha=0.5,edgecolor='dimgrey',zorder=2)
plt.scatter(4,dectrend_b2,color=cmap(0.7),s=70,alpha=1,edgecolor='dimgrey',zorder=2)

plt.scatter(5,dectrend_no1,color=cmap(0.9),s=70,alpha=0.5,edgecolor='dimgrey',zorder=2)
plt.scatter(5,dectrend_no2,color=cmap(0.9),s=70,alpha=1,edgecolor='dimgrey',zorder=2)

plt.ylim([-0.7,0.7])
plt.yticks(np.arange(-1,1.1,0.2),map(str,np.round(np.arange(-1,1.1,0.2),2)))
plt.ylabel(r'\textbf{(Near-)Surface Temperature Trends [$^{\circ}$C/decade]')

plt.text(5.84,dectrend_e1,r'1979-2019',fontsize=8,color='k',ha='center',va='center')
plt.text(5.84,dectrend_e2,r'2000-2019',fontsize=8,color='k',ha='center',va='center')

plt.xlim([0.5,5.5])
plt.xticks(np.arange(1,6,1),[r'\textbf{ERA5}',r'\textbf{NCEP2}',
                             r'\textbf{GISTEMPv4}',r'\textbf{BEST}',
                             r'\textbf{NOAAGlobalTempv5}'],size=10)
if sliceperiod == 'FMA':
    plt.title(r'\textbf{FEBRUARY--MARCH--APRIL}',fontsize=19,color='k')
elif sliceperiod == 'April':
    plt.title(r'\textbf{APRIL}',fontsize=19,color='k')

plt.savefig(directoryfigure + 'ObsDatasetTrends_%s.png' % sliceperiod,dpi=300)
