"""
Plot monthly trends in the Arctic for the PCHZ100

Author    : Zachary M. Labe
Date      : 28 June 022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 

monq = 'April'
yearmin = 1979
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)

### Read in observations
obs = np.genfromtxt(directorydata + 'T2M_BoxNA/ERA5_T2M_BoxOfInterest_TimeSeries-%s.txt' % monq,unpack=True)
obsz = sts.zscore(obs)

### Read in teleconnections
pchZ100d = np.genfromtxt(directorydata + 'PCHZ100/PCHZ100_%s_%s_%s-%s_detrended.txt' % ('Z100',monq,
                                                                yearmin,yearmax),
                                                                unpack=True)
years_pchZ100 = pchZ100d[0,:]
pchZ100 = sts.zscore(pchZ100d[1,:])
    
### Calculate trends
slope_pchZ100,intercept_pchZ100,r_pchZ100,p_pchZ100,se_pchZ100 = sts.linregress(years,pchZ100)

### Calculate correlations
corr_pchZ100,pcorr_pchZ100 = sts.pearsonr(obsz,pchZ100)
print(corr_pchZ100,pcorr_pchZ100)
sys.exit()
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

rects = plt.bar(years,pchZ100,align='center',zorder=2,color='crimson',label=r'PCH-Z100',clip_on=False)
plt.plot(years,obsz,linestyle='--',dashes=(1,0.3),color='k',linewidth=2,clip_on=False,label='ERA5-T2M')

l = plt.legend(shadow=False,fontsize=7,loc='upper center',
            fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,0.07),
            labelspacing=1,columnspacing=1,handletextpad=0.4)

plt.xticks(np.arange(1980,2021,5),map(str,np.arange(1980,2021,5)),fontsize=7)
plt.yticks(np.arange(-5,5.5,0.5),map(str,np.arange(-5,5.5,0.5)),fontsize=7)
plt.xlim([1979,2020])
plt.ylim([-3,3])

plt.ylabel(r'\textbf{Standardized Index}',fontsize=8,color='k')
plt.title(r'\textbf{Observational Indices For %s}' % monq,fontsize=12,color='k')
plt.savefig(directoryfigure + 'PCHZ100_Indices-Obs_%s' % monq,dpi=300)


