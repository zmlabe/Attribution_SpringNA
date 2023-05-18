"""
Plot rolling correlation between different datasets

Author    : Zachary M. Labe
Date      : 15 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import scipy.stats as sts
import pandas as pd

### Read in data files from server
indexName = 'WPPRECT'
indexVar = 'P'

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/T2M_BoxNA/' 
directorydata2 = '/work/Zachary.Labe/Data/ClimateIndices/%s/' % indexName

### Parameters
variq = 'T2M'
months = 'FMA'

### Read data
years_era,tas_era = np.genfromtxt(directorydata + 'T2M_BoxNA_T2M_%s_1979-2019_detrended.txt' % months,
                     unpack=True)
years_era,al_era = np.genfromtxt(directorydata2 + '%s_%s_%s_1979-2019_detrended.txt' % (indexName,indexVar,months),
                     unpack=True)

years_20cr,tas_20cr = np.genfromtxt(directorydata + 'T2M_BoxNA-20CRv3_T2M_%s_1921-2015_detrended.txt' % months,
                     unpack=True)
years_20cr,al_20cr = np.genfromtxt(directorydata2 + '%s-20CRv3_%s_%s_1921-2015_detrended.txt' % (indexName,indexVar,months),
                     unpack=True)

### Calculate rolling correlation
window = 10
tas_erapd = pd.Series(tas_era)
al_erapd = pd.Series(al_era)
corr_era = tas_erapd.rolling(window).corr(al_erapd).to_numpy()

tas_20crpd = pd.Series(tas_20cr)
al_20crpd = pd.Series(al_20cr)
corr_20cr = tas_20crpd.rolling(window).corr(al_20crpd).to_numpy()

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
ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.35,clip_on=False)

plt.axhline(y=0,color='dimgrey',linestyle='-',linewidth=2)

plt.plot(years_20cr[:-window+1],corr_20cr[window-1:],linestyle='--',label=r'\textbf{20CRv3}',linewidth=2,
         color='teal',dashes=(1,0.3))
plt.plot(years_era[:-window+1],corr_era[window-1:],linestyle='-',label=r'\textbf{ERA5}',linewidth=2,
         color='maroon')
         
leg = plt.legend(shadow=False,fontsize=10,loc='upper center',
           fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0.5,0.08),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
for line,text in zip(leg.get_lines(),leg.get_texts()):
    text.set_color(line.get_color())

plt.xticks(np.arange(1920,2040,10),np.arange(1920,2040,10))
plt.yticks(np.round(np.arange(-1,1.1,0.1),2),np.round(np.arange(-1,1.1,0.1),2))
plt.xlim([1920,2020])
plt.ylim([-1,1])

plt.ylabel(r'\textbf{Correlation Coefficient}',fontsize=11,
                      color='dimgrey')
plt.xlabel(r'\textbf{First Year In Window}',fontsize=11,color='dimgrey')
plt.title(r'\textbf{Rolling Correlation between %s and T2M for %s}' % (indexName,months),
                    color='k',fontsize=15)
plt.tight_layout()        

### Save figure
plt.savefig(directoryfigure+'%s_rollingCorrelations_%s-T2M.png' % (indexName,months),dpi=300)
