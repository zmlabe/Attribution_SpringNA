"""
Calculate KDE plots of anomalies for hierarchy of models

Author    : Zachary M. Labe
Date      : 15 August 2022
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

### Loop through indices
indicesdir = ['NINO34','WPPRECT','NECPSST','AL','PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
indicesvari = ['SST','PRECT','SST','SLP','Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500']

### Parameters
yearmin = 1979
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)
mon_index = 'FMA'
mon_t2m = 'FMA'
splitmin = -1
splitmax = 1

def compositePhase(index,data,splitmin,splitmax):
    pos = []
    neg = []
    for e in range(data.shape[0]):
        for yr in range(data.shape[1]):
            if index[e,yr] >= splitmax:
                pos.append(data[e,yr])
            elif index[e,yr] <= splitmin:
                neg.append(data[e,yr])
    return np.asarray(pos),np.asarray(neg)

def calcBandwidth(data):
    ### Bandwidth cross-validation (jakevdp.github.io/blog/2013/12/01/kernel-density-esimation/)
    grid = GridSearchCV(KernelDensity(),
                        {'bandwidth':np.linspace(0.1,2.1,50)},
                         cv = 10)
    grid.fit(data[:,None])
    bandwidth = grid.best_params_['bandwidth']
    print('Bandwidth is ---> %s!' % bandwidth)
    return bandwidth

def kde_sklearn(x,grid,bandwidth,**kwargs):
    """kerndel density estimation with scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth,**kwargs)
    kde_skl.fit(x[:,np.newaxis])
    
    ### score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(grid[:,np.newaxis])
    
    return np.exp(log_pdf)

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
    
    spear_nat_tempq = np.load(directorydata + 'T2M_BoxNA_AllModels/ClimateIndices_coupledSPEAR-NATURAL-T2M_BoxNA_%s_1979-2019_raw.npz' % mon_t2m)
    spear_nat_temp = spear_nat_tempq['spear'][0]
    
    ### Read in SPEAR-MED
    if iii < 3:
        spear_med_indexq = np.load(directorydata + 'coupledSPEAR/ClimateIndices-Tropical_coupledSPEAR_%s_1979-2019_trendIncluded.npz' % mon_index)
        indices_med = spear_med_indexq['indices']
        indicesvari_med = spear_med_indexq['indicesvari']
        spear_med_index = spear_med_indexq['spear'][iii]
    else:
        spear_med_indexq = np.load(directorydata + 'coupledSPEAR/ClimateIndices_coupledSPEAR_%s_1979-2019_trendIncluded.npz' % mon_index)
        indicesvari_med = spear_med_indexq['indicesvari']
        spear_med_index = spear_med_indexq['spear'][iii-3]
    
    spear_med_tempq = np.load(directorydata + 'T2M_BoxNA_AllModels/ClimateIndices_coupledSPEAR-T2M_BoxNA_%s_1979-2019_raw.npz' % mon_t2m)
    spear_med_temp = spear_med_tempq['spear'][0]
    
    ### Composite onto indices
    spear_nat_pos,spear_nat_neg = compositePhase(spear_nat_index,spear_nat_temp,splitmin,splitmax)
    spear_med_pos,spear_med_neg = compositePhase(spear_med_index,spear_med_temp,splitmin,splitmax)
    
    ### Calculate bandwidth for KDE
    spear_nat_pos_band = calcBandwidth(spear_nat_pos)
    spear_nat_neg_band = calcBandwidth(spear_nat_neg)
    spear_med_pos_band = calcBandwidth(spear_med_pos)
    spear_med_neg_band = calcBandwidth(spear_med_neg)
    
    ### Calculate KDE
    grid = np.arange(-15,15.1,0.1)
    spear_nat_pos_kde = kde_sklearn(spear_nat_pos,grid,spear_nat_pos_band)
    spear_nat_neg_kde = kde_sklearn(spear_nat_neg,grid,spear_nat_neg_band)
    
    spear_med_pos_kde = kde_sklearn(spear_med_pos,grid,spear_med_pos_band)
    spear_med_neg_kde = kde_sklearn(spear_med_neg,grid,spear_med_neg_band)
    
    ### Calculate KS test
    stat_nat,p_nat = sts.ks_2samp(spear_nat_pos,spear_nat_neg)
    stat_med,p_med = sts.ks_2samp(spear_med_pos,spear_med_neg)

    ###############################################################################
    ###############################################################################
    ###############################################################################    
    ### Plot PDFs
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})
    
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
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_color('dimgrey')
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
    
    plt.axvline(x=0,linestyle='-',linewidth=2,color='dimgrey')
      
    plt.plot(grid,spear_med_pos_kde,color='crimson',linewidth=2.1,
             linestyle='-',label=r'\textbf{SPEAR-MED [+%s]}' % indicesdir[iii])
    plt.plot(grid,spear_med_neg_kde,color='crimson',linewidth=1.4,
             linestyle='--',dashes=(1,0.5),label=r'\textbf{SPEAR-MED [-%s]}' % indicesdir[iii])  
      
    plt.plot(grid,spear_nat_pos_kde,color='deepskyblue',linewidth=2.1,
             linestyle='-',label=r'\textbf{SPEAR-MED-NATURAL [+%s]}' % indicesdir[iii])
    plt.plot(grid,spear_nat_neg_kde,color='deepskyblue',linewidth=1.4,
             linestyle='--',dashes=(1,0.5),label=r'\textbf{SPEAR-MED-NATURAL [-%s]}' % indicesdir[iii])
    
    plt.yticks(np.arange(0,1,0.1),list(map(str,np.round(np.arange(0,1,0.1),2))),
               fontsize=6)
    plt.xticks(np.arange(-15,15.1,0.5),list(map(str,np.arange(-15,15.1,0.5))),
               fontsize=6) 
    plt.xlim([-6,6])
    plt.ylim([0,0.4])
    
    l = plt.legend(shadow=False,fontsize=10,loc='upper center',
               fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0.5,1.17),
               labelspacing=0.2,columnspacing=1,handletextpad=0.4)
    for text in l.get_texts():
        text.set_color('k')
        
    plt.text(-6,0.35,r'\textit{p}=%s' % (format(np.round(p_nat,8),'f')),color='deepskyblue',fontsize=8)
    plt.text(-6,0.33,r'\textit{p}=%s' % (format(np.round(p_med,8),'f')),color='crimson',fontsize=8)
    
    plt.ylabel(r'\textbf{Density}',color='k',fontsize=8)  
    plt.xlabel(r'\textbf{T2M-BoxNA Anomalies [$^{\circ}$C]}',color='k',fontsize=8)
    
    plt.savefig(directoryfigure + 'KDEs_%s-Index-%s_T2M-%s_SPEAR.png' % (indicesdir[iii],mon_index,mon_t2m),dpi=300)
