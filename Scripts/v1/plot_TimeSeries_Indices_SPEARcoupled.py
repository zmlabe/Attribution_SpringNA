"""
Calculate timeseries for different climate indices

Author    : Zachary M. Labe
Date      : 16 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts
import cmocean
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/CoupledSPEAR/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/'
directorydata2 = '/work/Zachary.Labe/Data/ClimateIndices/'

### Loop through indices
indicesdir = ['NINO34','WPPRECT','NECPSST','AL','NAO','NPI','PCHT2M','SHI']
indicesvari = ['SST','P','SST','SLP','SLP','SLP','T2M','SLP']

### Parameters
yearmin = 1979
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)
mon_index = 'FMA'

for iii in range(len(indicesdir)):

    ### Read in SPEAR-NATURAL
    if iii < 3:
        spear_nat_indexq = np.load(directorydata + 'coupledSPEAR-NATURAL/ClimateIndices-Tropical_coupledSPEAR-NATURAL_%s_1979-2019.npz' % mon_index)
        indices_nat = spear_nat_indexq['indices']
        indicesvari_nat = spear_nat_indexq['indicesvari']
        spear_nat_index = spear_nat_indexq['spear'][iii]
    else:
        spear_nat_indexq = np.load(directorydata + 'coupledSPEAR-NATURAL/ClimateIndices_coupledSPEAR-NATURAL_%s_1979-2019.npz' % mon_index)
    indices_nat = spear_nat_indexq['indices']
    indicesvari_nat = spear_nat_indexq['indicesvari']
    spear_nat_index = spear_nat_indexq['spear'][iii-3]
    
    ### Read in SPEAR-MED
    if iii < 3:
        spear_med_indexq = np.load(directorydata + 'coupledSPEAR/ClimateIndices-Tropical_coupledSPEAR_%s_1979-2019_trendIncluded.npz' % mon_index)
        indices_med = spear_med_indexq['indices']
        indicesvari_med = spear_med_indexq['indicesvari']
        spear_med_index = spear_med_indexq['spear'][iii]
    else:
        spear_med_indexq = np.load(directorydata + 'coupledSPEAR/ClimateIndices_coupledSPEAR_%s_1979-2019_trendIncluded.npz' % mon_index)
        indicesvari_med = spear_med_indexq['indicesvari'][np.array([0,2,3,4,8])]
        spear_med_index = spear_med_indexq['spear'][np.array([0,2,3,4,8])][iii-3]
        
    ### Read in obs
    if indicesdir[iii] == 'NAO':
        indexd = np.genfromtxt(directorydata2 + '%s/%smodified_%s_%s_%s-%s.txt' % (indicesdir[iii],indicesdir[iii],'%s' % indicesvari[iii],mon_index,yearmin,yearmax),unpack=True)
    else:
        indexd = np.genfromtxt(directorydata2 + '%s/%s_%s_%s_%s-%s.txt' % (indicesdir[iii],indicesdir[iii],'%s' % indicesvari[iii],mon_index,yearmin,yearmax),unpack=True)
    years_index = indexd[0,:]
    indexo = sts.zscore(indexd[1,:])

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
    
    corrs = np.empty((len(spear_med_index)))
    for e in range(len(spear_med_index)):
        plt.plot(years,spear_med_index[e,:],linestyle='--',dashes=(1,0.3),color='darkgrey',linewidth=0.8,
                 alpha=0.4,clip_on=False)
        corrs[e] = sts.pearsonr(indexo,spear_med_index[e,:])[0]
    plt.plot(years,np.nanmean(spear_med_index,axis=0),linestyle='-',color='dimgrey',linewidth=2,clip_on=False,
             label=r'\textbf{SPEAR_MED}')
    plt.plot(years,np.nanmean(spear_nat_index,axis=0),linestyle='--',color='darkblue',linewidth=2,clip_on=False,
             label=r'\textbf{SPEAR_MED_NATURAL}',dashes=(1,0.3))
    plt.plot(years,indexo,linestyle='-',color='r',linewidth=2,clip_on=False,
             marker='o',label=r'\textbf{ERA5}')
    
    l = plt.legend(shadow=False,fontsize=7,loc='upper center',
                fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,0.07),
                labelspacing=1,columnspacing=1,handletextpad=0.4)
    
    plt.xticks(np.arange(1980,2021,5),map(str,np.arange(1980,2021,5)),fontsize=7)
    plt.yticks(np.arange(-5,5.5,0.5),map(str,np.arange(-5,5.5,0.5)),fontsize=7)
    plt.xlim([1979,2020])
    plt.ylim([-4,4])
    
    meancorr = np.round(np.mean(corrs),2)
    maxcorr = np.round(np.max(corrs),2)
    mincorr = np.round(np.min(corrs),2)
    plt.text(1979,-5.2,r'R=%s [%s,%s]' % (meancorr,mincorr,maxcorr),fontsize=10,color='k')
    
    plt.ylabel(r'\textbf{Standardized Index -- %s -- NOT DETRENDED}' % indicesdir[iii],fontsize=8,color='k')
    plt.title(r'\textbf{%s for %s}' % (indicesdir[iii],mon_index),fontsize=12,color='k')
    plt.savefig(directoryfigure + 'Index-%s_SPEAR_%s.png' % (indicesdir[iii],mon_index),dpi=300)
