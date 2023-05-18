"""
Plot box plots of trends of pattern correlations for different regions

Author    : Zachary M. Labe
Date      : 20 June 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import itertools

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/validationSPEAR/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/AMIPs/' 

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

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m"]
scenario = ['amip_1880s_rf','amip_obs_rf','spear']
model_1880s_rf = ['ECHAM5','ESRL-CAM5']
model_obs_rf = ['ECHAM5','ESRL-CAM5']
model_spear = ['spear']
experi = ['ECHAM5','ESRL-CAM5','SPEAR']
model = [model_1880s_rf,model_obs_rf,model_spear]
modelunravel = list(itertools.chain(*model))
variq = 'Z50'
slicemonth = 'April'
trendper = 'all'
slicenan = 'nan'
datareader = True
years = np.arange(1979,2019+1,1)

### Read in data
yrmin = 1979
yrmax = 2019
directorycorr = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/PatternCorrelations/'
corr = np.load(directorycorr + 'PatternCorrelations_Regional_%s_%s_validationSPEAR.npz' % (variq,slicemonth),
                allow_pickle=True)
echam = corr['echam']
cam = corr['cam']
spear = corr['spear']

### Look at different regions
regions = ['global','focusarea','arctic','NH','NHextra','NA','US']
regionsname = ['GLOBAL','REGION OF INTEREST','ARCTIC','NORTHERN HEMISPHERE','EXTRATROPICS -- NH','NORTH AMERICA','UNITED STATES']
plot6 = [0,1,2,3,5,6]

###############################################################################
###############################################################################
###############################################################################
### Graph for accuracy

fig = plt.figure()
for plo in range(len(plot6)):
    ax = plt.subplot(2,3,plo+1)
    
    plotdata = [echam[plo],cam[plo],spear[plo]]
    
    adjust_spines(ax, ['left', 'bottom'])
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_linewidth(2)
    ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
    ax.tick_params(axis="x",which="both",bottom = False,top=False,
                    labelbottom=False)
    
    ax.yaxis.grid(zorder=1,color='darkgrey',alpha=0.7,clip_on=False,linewidth=0.5)
    
    def set_box_color(bp, color):
        plt.setp(bp['boxes'],color=color)
        plt.setp(bp['whiskers'], color=color,linewidth=1.5)
        plt.setp(bp['caps'], color='w',alpha=0)
        plt.setp(bp['medians'], color='w',linewidth=1)
    
    positionsq = np.arange(len(experi))
    bpl = plt.boxplot(plotdata,positions=positionsq,widths=0.6,
                      patch_artist=True,sym='')
    
    # Modify boxes
    cp= 'maroon'
    set_box_color(bpl,cp)
    plt.plot([], c=cp, label=r'\textbf{PATTERN CORRELATION}',clip_on=False)
        
    for i in range(len(plotdata)):
        y = plotdata[i]
        x = np.random.normal(positionsq[i], 0.04, size=len(y))
        plt.plot(x, y,color='teal', alpha=0.5,zorder=10,marker='.',linewidth=0,markersize=5,markeredgewidth=0,clip_on=False)
     
    if any([plo==0,plo==3]):
        plt.yticks(np.arange(-1,1.1,0.2),list(map(str,np.round(np.arange(-1,1.1,0.2),2))),
                    fontsize=6) 
        plt.ylim([-1,1])
    else:
        plt.yticks(np.arange(-1,1.1,0.2),list(map(str,np.round(np.arange(-1,1.1,0.2),2))),
                    fontsize=6) 
        plt.ylim([-1,1])
        ax.axes.yaxis.set_ticklabels([])
    
    if plo>2:
        plt.text(-0.35,-1.1,r'\textbf{ECHAM5}',fontsize=5,color='dimgrey',
                  ha='left',va='center')
        plt.text(1.05,-1.1,r'\textbf{ESRL-CAM5}',fontsize=5,color='dimgrey',
                  ha='center',va='center')
        plt.text(2.3,-1.1,r'\textbf{SPEAR}',fontsize=5,color='dimgrey',
                  ha='right',va='center')
        
    plt.title(r'\textbf{%s}' % regionsname[plo],fontsize=8,color='dimgrey')
    plt.text(-0.85,1.05,r'\textbf{[%s]}' % letters[plo],color='k',fontsize=6)
    
    if any([plo==0,plo==3]):
        plt.ylabel(r'\textbf{Pattern Correlation -- %s}' % variq,color='k',fontsize=7)

plt.savefig(directoryfigure + 'patternCorrelations_regions_validationSPEAR_%s_%s_%s.png' % (variq,slicemonth,trendper),dpi=300)
