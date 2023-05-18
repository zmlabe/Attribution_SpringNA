"""
Calculate Aleutian Low Index

Author    : Zachary M. Labe
Date      : 26 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_ERA5_monthly as ERA
import calc_DetrendData as DT
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/Users/zlabe/Documents/Research/Attribution_SpringNA/Figures/' 
directorydata = '/Users/zlabe/Data/ERA5/' 

### Create months to loop through
monthsall = ['MAM']

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
variq = 'SLP'
sliceperiod = monthsall[0]
yearsall = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
yearmin = 1979
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)

### Read data
latobs,lonobs,varn = ERA.read_ERA5_monthly(variq,directorydata,sliceperiod,yearsall,
                                               sliceshape,addclimo,slicenan,False)
lon2,lat2 = np.meshgrid(lonobs,latobs)

### Read only 1979-2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]
var = varn[yearq,:,:]

### Calculate anomalies
anoms = np.asarray(calc_anomalies(years,var))

### Detrend data
vardt = DT.detrendDataR(anoms,'surface','monthly')
# vardt = anoms

### Calculate AL
la1 = 45
la2 = 65
lo1 = 160
lo2 = 200
latq = np.where((latobs >= la1) & (latobs <= la2))[0]
lonq = np.where((lonobs >= lo1) & (lonobs <= lo2))[0]
anomlon = vardt[:,:,lonq]
anoms = anomlon[:,latq,:]
lat2sq = lat2[latq,:]
lat2s = lat2sq[:,lonq]
meanz = UT.calc_weightedAve(anoms,lat2s)
mean = sts.zscore(meanz)
                                                             
print('\n========Calculated AL - %s =======\n' % monthsall[0])

slope, intercept, r, p, se = sts.linregress(years,mean)
trendline = slope*years + intercept

slopeRE, interceptRE, rRE, pRE, seRE = sts.linregress(years[:35],mean[:35])
trendlineRE = slopeRE*years[:35] + interceptRE

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
plt.plot(years[:35],trendlineRE,linewidth=2,color='steelblue',clip_on=False)

# plt.axhline(y=mean81,xmin=0.05,xmax=0.723,color='k',linewidth=2,
#             linestyle='--')
# plt.axhline(y=mean91,xmin=0.28,xmax=0.955,color='k',linewidth=2,
#             linestyle='--')

plt.xticks(np.arange(1980,2040,10),np.arange(1980,2040,10))
plt.yticks(np.round(np.arange(-9,9.1,0.5),2),np.round(np.arange(-9,9.1,0.5),2))
plt.xlim([1979,int(years[-1])])
plt.ylim([-2.5,2.5])

plt.ylabel(r'\textbf{SLP [hPa]}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{%s - ALEUTIAN LOW INDEX}' % sliceperiod,
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'%s_ALindex_1979-2019.png' % sliceperiod,dpi=300)

