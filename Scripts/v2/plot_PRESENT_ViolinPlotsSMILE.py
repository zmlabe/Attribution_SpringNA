"""
Plot box plots of trends over different time periods compared to SMILE for AGU

Author    : Zachary M. Labe
Date      : 8 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import scipy.stats as sts
import palettable.cubehelix as cm
import palettable.scientific.sequential as sss
import palettable.cartocolors.qualitative as cc
import cmocean as cmocean
import cmasher as cmr

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/Presentations/' 
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/MMLEA/'

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'FMA'

### Read in data
yrmin = 1979
yrmax = 2020
trend_obs = np.loadtxt(directoryoutput + 'Slopes_%s_obs_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_lens1 = np.loadtxt(directoryoutput + 'Slopes_%s_lens1_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_lens2 = np.loadtxt(directoryoutput + 'Slopes_%s_lens2_%s-%s.txt' % (sliceperiod,yrmin,yrmax))

trend_miroc = np.loadtxt(directoryoutput + 'Slopes_%s_MIRO6_LE_LOWS_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_smhi = np.loadtxt(directoryoutput + 'Slopes_%s_SMHI_LE_LOWS_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_mpi = np.loadtxt(directoryoutput + 'Slopes_%s_MPI_ESM12_LE_%s-%s.txt' % (sliceperiod,yrmin,yrmax))

trend_gfdlem = np.loadtxt(directoryoutput + 'Slopes_%s_gfdlem_%s-%s.txt' % (sliceperiod,yrmin,yrmax))
trend_gfdlcm = np.loadtxt(directoryoutput + 'Slopes_%s_gfdlc_%s-%s.txt' % (sliceperiod,yrmin,yrmax))

### Calculate decadal trend
decobs = trend_obs * 10
declens1= trend_lens1 * 10
declens2 = trend_lens2 * 10

decmiroc = trend_miroc * 10
decsmhi = trend_smhi * 10
decmpi = trend_mpi * 10

decgfdlem = trend_gfdlem * 10
decgfdlcm = trend_gfdlcm * 10

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

fig = plt.figure(figsize=(4,7))
axb = plt.subplot(111)

axb.spines['top'].set_color('none')
axb.spines['right'].set_color('none')
axb.spines['left'].set_color('none')
axb.spines['bottom'].set_color('none')
axb.tick_params('both',length=4,width=1,which='major',color='darkgrey')
axb.yaxis.grid(zorder=1,color='darkgrey',alpha=1,linewidth=1,clip_on=False)
plt.axvline(x=0,color='darkgrey',linestyle='-',linewidth=1,zorder=1)

datum = [decgfdlem,decgfdlcm,declens1,declens2,decmiroc,decsmhi,decmpi]

vp = plt.violinplot(datum,showmeans=False,showmedians=True,
                    vert=False,widths=0.7,showextrema=True)
plt.axvline(decobs,color='crimson',linestyle='--',dashes=(1,0.3),linewidth=2)
plt.setp(axb,yticks=[y+1 for y in range(len(datum))],
                     yticklabels=['GFDL_ESM2M_LE (30)','GFDL_CM3_LE (20)','CESM1_LE (40)','CESM2_LE (100)','MIROC6_LE (50)','SMHI_LE (50)', 'MPI_ESM1.2LR_LE (30)'])

positionsq = np.array(np.arange(1,8,1))
for i in range(len(datum)):
    y = datum[i]
    x = np.random.normal(positionsq[i], 0.02, size=len(y))
    plt.plot(y,x,color='w',alpha=0.5,zorder=10,marker='.',linewidth=0,markersize=5,markeredgewidth=0,clip_on=False)
                     
for i in vp['bodies']:
    i.set_edgecolor('w') 
    i.set_linewidth(0.4)
vp['cbars'].set_color('w')
vp['cmaxes'].set_color('w')
vp['cmins'].set_color('w')
vp['cmedians'].set_color('w')
vp['cmedians'].set_linewidth(2)
vp['cmaxes'].set_linewidth(2)        
vp['cmins'].set_linewidth(2)       
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')  

color = cmr.rainforest(np.linspace(0.00,0.8,len(datum)))
for i,c in zip(range(len(datum)),color):       
    vp['bodies'][i].set_facecolor(c)    
vp['bodies'][0].set_alpha(1)
vp['bodies'][1].set_alpha(1)  
vp['bodies'][2].set_alpha(1)   
vp['bodies'][3].set_alpha(1)
vp['bodies'][4].set_alpha(1)  
vp['bodies'][5].set_alpha(1)
vp['bodies'][6].set_alpha(1) 

plt.xticks(np.arange(-5,5.01,0.5),map(str,np.arange(-5,5.01,0.5)))
plt.yticks(color='w')
plt.xlim([-0.5,1.2])  
plt.xlabel(r'\textbf{T2M Trends [$^{\circ}$C/Decade] for April %s-%s}' % (yrmin,yrmax),color='w',size=7)

axb.xaxis.set_ticks_position('bottom')
axb.yaxis.set_ticks_position('left')

plt.text(decobs,7.8,r'\textbf{ERA5}',size=11,color='crimson',ha='center',va='center')

plt.tight_layout()
plt.savefig(directoryfigure + 'PRESENT_TrendViolins_CESM_%s_%s-%s.png' % (sliceperiod,yrmin,yrmax),dpi=300)







