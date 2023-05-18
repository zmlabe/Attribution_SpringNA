"""
Plot histogram of indices for correlations between index and obs in AMIPS
for the tropics

Author    : Zachary M. Labe
Date      : 29 November 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/AMIPs/PDFs/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/'
directorydata3 = '/work/Zachary.Labe/Data/ClimateIndices/' 

### Loop through indices
indicesdir = ['NINO34','WPPRECT','NECPSST']
variqn = ['SST','PRECT','SST']
models = ['spear']
modelname = ['spear']
scenario = ['spear_obs_rf']

### Parameters
yearmin = 1979
yearmax = 2020
years = np.arange(yearmin,yearmax+1,1)
mon_index = 'April'
mon_t2m = 'April'

### Coldest AMIP ensemble member
loadens_amip,min_amip = np.load(directorydata2 + 'EnsembleMembers_SimilarTrends_AMIPSPEAR_T2M_%s_1979-2020.npy' % mon_t2m,
                        allow_pickle=True)

for iii in range(len(indicesdir)):
    
    ### Read in AMIP data
    amip = np.load(directorydata2 + 'ClimateIndices-Tropical_AMIPSPEAR_%s_1979-2020.npz' % mon_index)
    indices = amip['indices']
    indicesvari = amip['indicesvari']
    yearsamip = amip['years']
    
    spearindex = amip['spear'][iii]
    
    ### Read in AMIP T2M data
    amiptas = np.load(directorydata2 + 'ClimateIndices-T2M_BoxNA_AMIPSPEAR_%s_1979-2020.npz' % mon_t2m)
    indicestas = amiptas['indices']
    indicesvaritas = amiptas['indicesvari']
    yearsamiptas = amiptas['years']
    
    speartas = amiptas['spear'][0]
    
    ### Calculate correlation
    corrall = np.empty((len(speartas)))
    pall = np.empty((len(speartas)))
    for ens in range(len(speartas)):
        corrall[ens],pall[ens] = sts.pearsonr(speartas[ens,:],spearindex[ens,:])
        
    corrmin,pmin = sts.pearsonr(speartas[min_amip,:],spearindex[min_amip,:])
        
    ###########################################################################
    ### Observations
    obs = np.genfromtxt(directorydata3 + 'T2M_BoxNA/T2M_BoxNA_%s_%s_%s-%s_detrended.txt' % ('T2M',mon_t2m,
                                                                    yearmin,yearmax),
                                                                    unpack=True)
    obsz = sts.zscore(obs[1,:])
    
    ### Read in indices
    indexallname = ['NINO34','WPPRECT','NECPSST']
    variqn = ['SST','P','SST']
    

    if any([indexallname[iii]=='PNA',indexallname[iii]=='NAO']):
        indexq = np.genfromtxt(directorydata3 + '%s/%smodified_%s_%s_%s-%s_detrended.txt' % (indexallname[iii],
                                                                        indexallname[iii],variqn[iii],mon_index,
                                                                        yearmin,yearmax),
                                                                        unpack=True)
    else:
        indexq = np.genfromtxt(directorydata3 + '%s/%s_%s_%s_%s-%s_detrended.txt' % (indexallname[iii],
                                                                        indexallname[iii],variqn[iii],mon_index,
                                                                        yearmin,yearmax),
                                                                        unpack=True)
    indexall = sts.zscore(indexq[1,:])
    print('Read in data for ---> %s' % indexallname[iii])
        
    ### Correlation observations
    corrobs,pobs = sts.pearsonr(obsz,indexall)

    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Create plot for histograms of slopes
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
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
    adjust_spines(ax, ['left','bottom'])            
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none') 
    ax.spines['bottom'].set_color('dimgrey')
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2) 
    ax.tick_params('both',length=5.5,width=2,which='major',color='dimgrey')  
    ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.35)
    
    ### Plot histograms      
    weights_corrall = np.ones_like(corrall)/len(corrall)
    n_corrall, bins_corrall, patches_corrall = plt.hist(corrall,bins=np.arange(-1,1.01,0.1),
                                            density=False,alpha=0.5,
                                            label=r'\textbf{%s}' % indicesdir[iii],
                                            weights=weights_corrall,zorder=5)
    for i in range(len(patches_corrall)):
        patches_corrall[i].set_facecolor('teal')
        patches_corrall[i].set_edgecolor('white')
        patches_corrall[i].set_linewidth(0.5)
        
    leg = plt.legend(shadow=False,fontsize=7,loc='upper center',
            bbox_to_anchor=(0.5,1.1),fancybox=True,ncol=3,frameon=False,
            handlelength=3,handletextpad=1)
    
    plt.axvline(x=0,linewidth=2,color='dimgrey',linestyle='-',zorder=11)
    plt.axvline(x=corrobs,linewidth=4,color='maroon',linestyle='--',dashes=(1,0.3),zorder=12)
    plt.axvline(x=corrmin,linewidth=2,color='gold',linestyle='--',dashes=(1,0.3),zorder=12)
    
    plt.ylabel(r'\textbf{PROPORTION}',fontsize=10,color='k')
    plt.xlabel(r'\textbf{R} [AMIP-SPEAR - %s -- 1979-2020]' % mon_index,fontsize=10,color='k')
    plt.yticks(np.arange(0,1.1,0.1),map(str,np.round(np.arange(0,1.1,0.1),2)),size=6)
    plt.xticks(np.arange(-1,1.1,0.1),map(str,np.round(np.arange(-1,1.1,0.1),2)),size=6)
    plt.xlim([-1,1])   
    plt.ylim([0,0.4])
        
    plt.savefig(directoryfigure + 'Histogram_r_1979-2020_AMIPSPEAR_%s_%s.png' % (indicesdir[iii],mon_index),
                dpi=300)
    
    
