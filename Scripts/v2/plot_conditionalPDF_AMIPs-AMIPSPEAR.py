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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/AMIPs/PDFs/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/'

### Loop through indices
indicesdir = ['NPOcustomZ','ABI','NPO','AL','PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
models = ['spear']
modelname = ['spear']
scenario = ['spear_obs_rf']

### Parameters
yearmin = 1979
yearmax = 2020
years = np.arange(yearmin,yearmax+1,1)
mon_index = 'FMA'
mon_t2m = 'FMA'

for iii in range(len(indicesdir)):
    
    ### Read in AMIP data
    amip = np.load(directorydata2 + 'ClimateIndices_AMIPSPEAR_%s_1979-2020.npz' % mon_index)
    indices = amip['indices']
    indicesvari = amip['indicesvari']
    yearsamip = amip['years']
    
    spear = amip['spear'][iii]
    
    ### Read in AMIP T2M data
    amiptas = np.load(directorydata2 + 'ClimateIndices-T2M_BoxNA_AMIPSPEAR_%s_1979-2020.npz' % mon_t2m)
    indicestas = amiptas['indices']
    indicesvaritas = amiptas['indicesvari']
    yearsamiptas = amiptas['years']
    
    speartas = amiptas['spear'][0]
    
    def compositePhase(data,index,splitmin,splitmax):
        pos = []
        neg = []
        for e in range(data.shape[0]):
            for yr in range(data.shape[1]):
                if index[e,yr] > splitmax:
                    pos.append(data[e,yr])
                elif index[e,yr] < splitmin:
                    neg.append(data[e,yr])
        return pos,neg
    
    ### Composite onto indices
    splitmin = -0.5
    splitmax = 0.5
    spear_pos,spear_neg = compositePhase(speartas,spear,splitmin,splitmax)
    
    def calcPDFs(pos,neg):
        num_bins = np.arange(-4,4.01,0.01)
        
        m_pos,s_pos = sts.norm.fit(pos)
        pdf_pos = sts.norm.pdf(num_bins,m_pos,s_pos)
        
        m_neg,s_neg = sts.norm.fit(neg)
        pdf_neg = sts.norm.pdf(num_bins,m_neg,s_neg)
        
        return pdf_pos,pdf_neg,num_bins
                    
    ### Calculate PDFs
    spear_pdf_pos,spear_pdf_neg,bins  = calcPDFs(spear_pos,spear_neg)
    
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
    
    plt.plot(bins,spear_pdf_pos,color='goldenrod',linewidth=2.1,
             linestyle='-',label=r'\textbf{SPEAR-MED-AMIP [+%s]}' % indices[iii],clip_on=False)
    plt.plot(bins,spear_pdf_neg,color='goldenrod',linewidth=1.4,
             linestyle='--',dashes=(1,0.5),label=r'\textbf{SPEAR-MED-AMIP [-%s]}' % indices[iii],clip_on=False)
    
    plt.yticks(np.arange(0,1,0.1),list(map(str,np.round(np.arange(0,1,0.1),2))),
               fontsize=6)
    plt.xticks(np.arange(-5,5.1,0.5),list(map(str,np.arange(-5,5.1,0.5))),
               fontsize=6) 
    plt.xlim([-4,4])
    plt.ylim([0,0.5])
    
    l = plt.legend(shadow=False,fontsize=10,loc='upper center',
               fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,1.17),
               labelspacing=0.2,columnspacing=1,handletextpad=0.4)
    for text in l.get_texts():
        text.set_color('k')
    
    plt.ylabel(r'\textbf{Density}',color='k',fontsize=8)  
    plt.xlabel(r'\textbf{T2M-BoxNA [$^{\circ}$C]}',color='k',fontsize=8)
    
    plt.savefig(directoryfigure + 'PDFs_%s-Index-%s_T2M-%s_AMIPs.png' % (indices[iii],mon_index,mon_t2m),dpi=300)
