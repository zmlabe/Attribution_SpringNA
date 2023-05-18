"""
Calculate regression maps for loops of indices

Author    : Zachary M. Labe
Date      : 17 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
import calc_DetrendData as DT
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import cmocean

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/ObsIndicesAMIPS/' 
directorydataObs = '/work/Zachary.Labe/Data/ClimateIndices/'
directorydataObsT2Mbox = '/work/Zachary.Labe/Data/ClimateIndices/T2M_BoxNA/'

### Parameters
yearsall = np.arange(1979,2021+1,1)
yearmin = 1979
yearmax = 2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

### Loop through indices
indicesdir = ['PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI','NINO34','WPPRECT','NECPSST','SND_BoxNA']
indices = ['PNAmodified','NAOmodified','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI','NINO34','WPPRECT','NECPSST','SND_BoxNA']
indicesvari = ['Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500','SST','P','SST','SND']
monthsloop = ['January','February','March','April']

### Read in indices from observations
robs = np.empty((len(indices),len(monthsloop),len(monthsloop)))
pobs = np.empty((len(indices),len(monthsloop),len(monthsloop)))
for iii in range(len(indices)):
    for mmm in range(len(monthsloop)):
        for ooo in range(len(monthsloop)):
            index_mon_obs = monthsloop[mmm]
            t2m_mon_obs = monthsloop[ooo]
            
            ### Read in T2M box
            t2m_box_obs = np.genfromtxt(directorydataObsT2Mbox + 'ERA5_T2M_BoxOfInterest_TimeSeries-%s.txt' % t2m_mon_obs,unpack=True)
            yearclimoq = np.where((years>=1981) & (years<=2010))[0]
            years_BOXq = t2m_box_obs[yearclimoq]
            meanclimoBOX = np.nanmean(years_BOXq)
            stdclimoBOX = np.nanstd(years_BOXq)
            BOX = (t2m_box_obs - meanclimoBOX)/stdclimoBOX
            INDEXd = np.genfromtxt(directorydataObs + '%s/%s_%s_%s_%s-%s_detrended.txt' % (indicesdir[iii],indices[iii],indicesvari[iii],index_mon_obs,yearmin,yearmax),unpack=True)
            years_INDEX = INDEXd[1,:]
            yearclimoq = np.where((years>=1981) & (years<=2010))[0]
            years_INDEXq = years_INDEX[yearclimoq]
            meanclimo = np.nanmean(years_INDEXq)
            stdclimo = np.nanstd(years_INDEXq)
            INDEX = (years_INDEX - meanclimo)/stdclimo
            
            ### Calculate correlations
            robs[iii,mmm,ooo],pobs[iii,mmm,ooo] = sts.pearsonr(BOX,INDEX)
            
###############################################################################
###############################################################################
###############################################################################
for iii in range(len(indices)):
    fig = plt.figure()
    ax = plt.subplot(111)
    
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
                color='w')
    
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom=True)
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft='on')
    
    cs = plt.pcolormesh(robs[iii],shading='faceted',edgecolor='w',
                        linewidth=0.3,vmin=-1,vmax=1)
    
    for i in range(robs.shape[1]):
        for j in range(robs.shape[2]):
            if pobs[iii,i,j] <= 0.05:
                plt.text(j+0.5,i+0.5,r'\textbf{%+1.2f}' % robs[iii,i,j],fontsize=7,
                      color='k',va='center',ha='center')
    
    cs.set_cmap(cmocean.cm.balance)
    
    plt.yticks(np.arange(0.5,len(monthsloop)+0.5,1),monthsloop,ha='right',color='dimgrey',
                va='center')
    yax = ax.get_yaxis()
    yax.set_tick_params(pad=-2)
    plt.xticks(np.arange(0.5,len(monthsloop)+0.5,1),monthsloop,ha='center',color='dimgrey',
                va='center')
    xax = ax.get_xaxis()
    xax.set_tick_params(pad=8)
    plt.xlim([0,4])
    plt.ylim([0,4])
    
    plt.xlabel(r'\textbf{T2M Box}',color='dimgrey',fontsize=10)
    plt.ylabel(r'\textbf{%s}' % indicesdir[iii],color='dimgrey',fontsize=10)
    
    cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
    ticks = np.arange(-1,1.1,0.2)
    labels = list(map(str,np.round(np.arange(-1,1.1,0.2),2)))
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    cbar.ax.tick_params(axis='x', size=.001)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.set_label(r'\textbf{Correlation Coefficient}',
                    color='dimgrey',labelpad=3,fontsize=12)
    
    plt.tight_layout()
    plt.savefig(directoryfigure + 'HeatMap_LeadLag_ClimateIndices-Obs_%s' % indicesdir[iii],dpi=900)
