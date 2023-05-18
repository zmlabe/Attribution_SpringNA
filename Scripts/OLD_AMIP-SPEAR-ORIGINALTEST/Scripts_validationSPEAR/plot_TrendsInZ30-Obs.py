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
directorydata = '/Users/zlabe/Documents/Research/Attribution_SpringNA/Data/Indices/ClimateIndices/PCHZ30/' 

### Create months to loop through
monthsall = ['April']

years,meanz = np.genfromtxt(directorydata + 'PCHZ30_Z30_%s_1979-2019.txt' % monthsall[0],
                     unpack=True)
mean = sts.zscore(meanz)

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
plt.ylim([-3,3])

plt.ylabel(r'\textbf{Standardized}',fontsize=11,
                      color='dimgrey')
plt.title(r'\textbf{%s - PCH-Z30 INDEX}' % monthsall[0],
                    color='k',fontsize=17)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'%s_PCHZ30_1979-2019.png' % monthsall[0],dpi=300)

