"""
Plot timeseries of ABI

Author    : Zachary M. Labe
Date      : 16 November 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts
import statsmodels.api as sm

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 

monq = 'April'
yearmin = 1979
yearmax = 2020
years = np.arange(yearmin,yearmax+1,1)

### Read in observations
obs = np.genfromtxt(directorydata + 'T2M_BoxNA/T2M_BoxNA_%s_%s_%s-%s_detrended.txt' % ('T2M',monq,
                                                                yearmin,yearmax),
                                                                unpack=True)
obsz = sts.zscore(obs[1,:])

### Read in teleconnections
ABId = np.genfromtxt(directorydata + 'ABI/ABI_%s_%s_%s-%s_detrended.txt' % ('Z500',monq,
                                                                yearmin,yearmax),
                                                                unpack=True)
years_ABI = ABId[0,:]
ABI = sts.zscore(ABId[1,:])
    
### Calculate trends
slope_ABI,intercept_ABI,r_ABI,p_ABI,se_ABI = sts.linregress(years,ABI)

### Calculate correlations
corr_ABI,pcorr_ABI = sts.pearsonr(obsz,ABI)
print(corr_ABI,pcorr_ABI)

### Smooth NPO index
smoothed = sm.nonparametric.lowess(ABI,np.arange(years.shape[0]),
                                   frac=10/years.shape[0])

###############################################################################
###############################################################################
###############################################################################
### Plot monthly indices
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
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',
               labelbottom='off',bottom='off')

rects = plt.bar(years,ABI,align='center',zorder=2,color='crimson',label=r'ABI',clip_on=False)
plt.plot(years,obsz,linestyle='--',dashes=(1,0.3),color='k',linewidth=2,clip_on=False,label='ERA5-T2M')
plt.plot(years,smoothed[:,1],linestyle='-',linewidth=2,color='darkred',clip_on=False)

l = plt.legend(shadow=False,fontsize=7,loc='upper center',
            fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,0.07),
            labelspacing=1,columnspacing=1,handletextpad=0.4)

plt.xticks(np.arange(1980,2021,5),map(str,np.arange(1980,2021,5)),fontsize=7)
plt.yticks(np.arange(-5,5.5,0.5),map(str,np.arange(-5,5.5,0.5)),fontsize=7)
plt.xlim([1979,2020])
plt.ylim([-3,3])

plt.ylabel(r'\textbf{Standardized Index}',fontsize=8,color='k')
plt.title(r'\textbf{Observational Indices For %s}' % monq,fontsize=12,color='k')
plt.savefig(directoryfigure + 'ABI_Indices-Obs_%s_detrended' % monq,dpi=300)


