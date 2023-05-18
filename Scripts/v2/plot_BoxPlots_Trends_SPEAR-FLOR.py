"""
Plot box plots of trends over different time periods compared to all SPEAR/FLOR

Author    : Zachary M. Labe
Date      : 18 October 2022
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
sliceperiod = 'FMA'

### Read in data
yrmin = 1979
yrmax = 2020
trend_obs = np.loadtxt(directoryoutput + 'Slopes_%s_obs_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_spear = np.loadtxt(directoryoutput + 'Slopes_%s_spear_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_flor = np.loadtxt(directoryoutput + 'Slopes_%s_FLOR_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_spearnoaer = np.loadtxt(directoryoutput + 'Slopes_%s_noaerspear_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_spearnatural = np.loadtxt(directoryoutput + 'Slopes_%s_natural_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_amip = np.load(directorydataamip + 'AMIP_onlySPEAR_Trends_%s_T2M.npz' % sliceperiod,allow_pickle=True)['trend'][0][0][0]

### Calculate decadal trend
decobs = trend_obs * 10
decspear = trend_spear * 10
decflor= trend_flor * 10
decspearnoaer = trend_spearnoaer * 10
decspearnatural = trend_spearnatural * 10
decamip = trend_amip * 10

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
axb = plt.subplot(111)

axb.spines['top'].set_color('none')
axb.spines['right'].set_color('none')
axb.spines['left'].set_color('none')
axb.spines['bottom'].set_color('none')
axb.tick_params('both',length=4,width=1,which='major',color='darkgrey')
axb.yaxis.grid(zorder=1,color='darkgrey',alpha=1,linewidth=1,clip_on=False)
plt.axvline(x=0,color='darkgrey',linestyle='-',linewidth=1,zorder=1)

datum = [decamip,decspear,decflor,decspearnoaer,decspearnatural]

vp = plt.violinplot(datum,showmeans=False,showmedians=True,
                    vert=False,widths=0.6,showextrema=True)
plt.axvline(decobs,color='k',linestyle='--',dashes=(1,0.3),linewidth=2)
plt.setp(axb,yticks=[y+1 for y in range(len(datum))],
                     yticklabels=['SPEAR_MED_AMIP (30)','SPEAR_MED (30)','FLOR (30)','SPEAR_MED_NOAER (12)','SPEAR_MED_NATURAL (30)'])

positionsq = np.array(np.arange(1,6,1))
for i in range(len(datum)):
    y = datum[i]
    x = np.random.normal(positionsq[i], 0.02, size=len(y))
    plt.plot(y,x,color='k',alpha=0.5,zorder=10,marker='.',linewidth=0,markersize=10,markeredgewidth=0,clip_on=False)
                     
for i in vp['bodies']:
    i.set_edgecolor('dimgrey')  
vp['cbars'].set_color('k')
vp['cmaxes'].set_color('k')
vp['cmins'].set_color('k')
vp['cmedians'].set_color('k')
vp['cmedians'].set_linewidth(3)
vp['cmaxes'].set_linewidth(3)        
vp['cmins'].set_linewidth(3)       
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')          
vp['bodies'][0].set_facecolor('salmon')
vp['bodies'][1].set_facecolor('crimson')   
vp['bodies'][2].set_facecolor('goldenrod') 
vp['bodies'][3].set_facecolor('deepskyblue') 
vp['bodies'][4].set_facecolor('teal')      
vp['bodies'][0].set_alpha(0.8)
vp['bodies'][1].set_alpha(0.8)  
vp['bodies'][2].set_alpha(0.8) 
vp['bodies'][3].set_alpha(0.8)  
vp['bodies'][4].set_alpha(0.8)   

plt.xticks(np.arange(-5,5.01,0.5),map(str,np.arange(-5,5.01,0.5)))
plt.xlim([-1,1.5])  
plt.xlabel(r'\textbf{%s for %s-%s of T2M Trends [$^{\circ}$C/Decade]' % (sliceperiod,yrmin,yrmax),color='k',size=14)

axb.xaxis.set_ticks_position('bottom')
axb.yaxis.set_ticks_position('left')

plt.text(decobs,5.7,r'ERA5',size=11,color='k',ha='center',va='center')

plt.tight_layout()
plt.savefig(directoryfigure + 'TrendViolins_SPEAR-FLOR_%s_%s-%s.png' % (sliceperiod,yrmin,yrmax),dpi=300)

### Save ensemble members
spearAMIPens = decamip
spearens = decspear

uncertainty = 0.1
spear_q = np.where((spearAMIPens <= decobs+uncertainty))[0]
spearCoupled_q = np.where((spearens <= decobs+uncertainty))[0] 
SPEARargmin = np.argmin(spearAMIPens)
SPEARcoupledargmin = np.argmin(spearens)
np.save(directorydataamip + 'EnsembleMembers_SimilarTrends_AMIPSPEAR_%s_%s_%s-%s.npy' % (variq,sliceperiod,yrmin,yrmax),
         [spear_q,SPEARargmin])
np.save(directorydataamip + 'EnsembleMembers_SimilarTrends_CoupledSPEAR_%s_%s_%s-%s.npy' % (variq,sliceperiod,yrmin,yrmax),
         [spearCoupled_q,SPEARcoupledargmin])







