"""
Plot differences of climate baselines 1981-2010 and 1991-2020 by month
Author    : Zachary M. Labe
Date      : 16 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/' 

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'FMA'
years = np.arange(1979,2021+1,1)
sliceshape = 4
slicenan = 'nan'
addclimo = True
newon = True

### Read data
lat1,lon1,lev1,data = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                years,sliceshape,addclimo,
                                slicenan,'surface')

data = data[:-2,:,:]
years = years[:-2]

### Baseline - 1981-2010
yearqold = np.where((years >= 1981) & (years <= 2010))[0]
climold = np.nanmean(data[yearqold,:,:],axis=0)
anommean = data - climold

### Calculate region
la1 = 43
la2 = 60
lo1 = 240
lo2 = 295
lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]

meanlat = anommean[:,lat1q,:]
meanbox = meanlat[:,:,lon1q]
lon1a = lon1[lon1q]
lat1a = lat1[lat1q]
lon2q,lat2q = np.meshgrid(lon1a,lat1a)

### Calculate timeseries
mean = UT.calc_weightedAve(meanbox,lat2q)
# np.savetxt(directorydata2 + 'ERA5_T2M_BoxOfInterest_TimeSeries-January.txt',mean)
# sys.exit()
### Calculate statistics
slope, intercept, r, p, se = sts.linregress(years,mean)
trendline = slope*years + intercept

slopeRE, interceptRE, rRE, pRE, seRE = sts.linregress(years[-20:],mean[-20:])
trendlineRE = slopeRE*years[-20:] + interceptRE

### Decadal trend
dectrend = slope*10
print('\n\nDecadal trend is ------> %s!' % dectrend)

### Calculate 1981-2010 mean line
mean81 = np.nanmean(mean[2:2+30])
mean91 = np.nanmean(mean[12:12+30])

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
ax.tick_params(axis='x',labelsize=6,pad=1.5)
ax.tick_params(axis='y',labelsize=6,pad=1.5)

plt.plot(years,mean,linewidth=2.3,color='maroon',alpha=1,clip_on=False,
         marker='o',markersize=6)
plt.plot(years,trendline,linewidth=1,color='teal',clip_on=False)
plt.plot(years[-20:],trendlineRE,linewidth=2,color='steelblue',clip_on=False)

plt.axhline(y=mean81,xmin=0.05,xmax=0.723,color='k',linewidth=2,
            linestyle='--')
plt.axhline(y=mean91,xmin=0.28,xmax=0.955,color='k',linewidth=2,
            linestyle='--')

plt.xticks(np.arange(1980,2040,10),np.arange(1980,2040,10))
plt.yticks(np.round(np.arange(-5,5.1,0.5),2),np.round(np.arange(-5,5.1,0.5),2))
plt.xlim([1979,int(years[-1])])
plt.ylim([-5,5])

plt.ylabel(r'\textbf{$\bf{^\circ}$C}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{FEBRUARY-APRIL TEMPERATURE ANOMALIES}',
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'FMA_meanT_1979-2019.png',dpi=300)
