"""
Plot box plots of trends over different time periods compared to all SPEAR/FLOR
for different single forcing contributions

Author    : Zachary M. Labe
Date      : 20 October 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/TrendDistributions/' 
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/MMLEA/'
directorydataamip = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/'

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'April'

### Read in data
yrmin = 1979
yrmax = 2020
trend_obs = np.loadtxt(directoryoutput + 'Slopes_%s_obs_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_spear = np.loadtxt(directoryoutput + 'Slopes_%s_spear_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_spearnoaer = np.loadtxt(directoryoutput + 'Slopes_%s_noaerspear_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_spearnatural = np.loadtxt(directoryoutput + 'Slopes_%s_natural_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_amip = np.load(directorydataamip + 'AMIP_onlySPEAR_Trends_%s_T2M.npz' % sliceperiod,allow_pickle=True)['trend'][0][0][0]

### Calculate decadal trend
decobs = trend_obs * 10
decspear = trend_spear * 10
decspearnoaer = trend_spearnoaer * 10
decspearnatural = trend_spearnatural * 10
decamip = trend_amip * 10

### Calculate contributions
COP = np.nanmean(decamip) - np.nanmean(decspear)
GHG = np.nanmean(decspearnoaer) - np.nanmean(decspearnatural)
AER = np.nanmean(decspearnoaer) - np.nanmean(decspear)

ADD = GHG + AER
ANT = np.nanmean(decspear) - np.nanmean(decspearnatural)
linear = ANT - ADD

### Begin the plotting routine
bars = [GHG,AER,ADD,ANT,COP]

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

fig = plt.figure(figsize=(7,4))
ax = plt.subplot(111)
ax.yaxis.grid(zorder=1,color='darkgrey',alpha=1,clip_on=False)
adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(0)
ax.spines['left'].set_linewidth(2)
ax.spines['left'].set_color('dimgrey')
ax.tick_params(axis='both', direction='out',length=5.5,width=2,
               which='major',pad=3,color='dimgrey')

recdiff_masked = np.ma.masked_less_equal(bars,0)

rects1 = plt.bar(np.arange(len(bars))+0.5,bars,color='deepskyblue',
                edgecolor='deepskyblue',zorder=2)
rects2 = plt.bar(np.arange(len(bars))+0.5,recdiff_masked,color='crimson',
                edgecolor='crimson',zorder=2)

plt.axhline(decobs,color='darkgrey',linestyle='--',dashes=(1,0.3),linewidth=2)

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        if rect.get_height() > 0:
            ax.text(rect.get_x() + rect.get_width()/2., rect.get_height()+0.025,
                    r'\textbf{+%2.2f}' % height,ha='center',va='top',
                    zorder=11,color='dimgrey',fontsize=7)
        elif rect.get_height() <= 0:
            ax.text(rect.get_x() + rect.get_width()/2., rect.get_height()-0.01,
                    r'\textbf{%2.2f}' % height,ha='center',va='top',
                    zorder=11,color='dimgrey',fontsize=7)
autolabel(rects1)
autolabel(rects2)

xlabels = [r'\textbf{GHG-Forced}',r'\textbf{AER-Forced}',r'\textbf{GHG+AER}',r'\textbf{Anthropogenic}',
           r'\textbf{Prescribed SST}']
plt.xticks(np.arange(0,5,1)+0.5,xlabels,rotation=0,fontsize=6,color='k')
ylabels = map(str,np.round(np.arange(-2,2.1,0.2),2))
plt.yticks(np.arange(-2,2.1,0.2),ylabels,fontsize=6,color='k')
plt.ylim([-0.4,0.4])
plt.xlim([0,5])

plt.text(1.9,-0.14,r'\textbf{ERA5 Observed Trend}',color='darkgrey',fontsize=11)

plt.title(r'\textbf{FORCING CONTRIBUTIONS}',fontsize=23,color='w',alpha=1) 
plt.ylabel(r'\textbf{%s for %s-%s of T2M Trends [$^{\circ}$C/Decade]' % (sliceperiod,yrmin,yrmax),color='k',size=6)

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(directoryfigure + 'ForcingContributions_%s_%s-%s.png' % (sliceperiod,yrmin,yrmax),dpi=300)







