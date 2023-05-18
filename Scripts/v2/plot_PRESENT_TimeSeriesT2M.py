"""
Plot each teleconnection timeseries from AMIP-SPEAR for AGU

Author    : Zachary M. Labe
Date      : 8 December 2022
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
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/AGU_Poster/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/' 

monq = 'April'
yearsall = np.arange(1979,2021+1,1)
yearmin = 1979
yearmax = 2020
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

indicesdir = ['T2M_BoxNA']
indices = ['T2M_BoxNA']
indicesvari = ['T2M']
models = ['SPEAR']
modelname = ['spear']

### Read in new indices
for i in range(len(indices)):
    indexd = np.genfromtxt(directorydata + '%s/%s_%s_%s_%s-%s.txt' % (indicesdir[i],indices[i],'%s' % indicesvari[i],monq,yearmin,yearmax),unpack=True)
    years_index = indexd[0,:]
    indexo = sts.zscore(indexd[1,:])
    
    for m in range(len(models)):
        
        ### Read in AMIP data
        directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/AMIPs/'
        amip = np.load(directoryoutput + 'ClimateIndices-T2M_BoxNA_AMIPSPEAR_%s_1979-2020_trendIncluded.npz' % monq)
        amipi = amip['%s' % modelname[m]][i,:]

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
        
        fig = plt.figure(figsize=(8,3))
        ax = plt.subplot(111)
        adjust_spines(ax, ['left', 'bottom'])
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['bottom'].set_color('darkgrey')
        ax.spines['left'].set_linewidth(0)
        ax.spines['bottom'].set_linewidth(2)
        ax.tick_params('x',length=4,width=2,which='major',color='darkgrey',
                       labelbottom='off',bottom='off')
        ax.tick_params('y',length=0,width=0,which='major',color='k',
                       labelbottom='off',bottom='off')
        
        corrs = np.empty((len(amipi)))
        for e in range(len(amipi)):
            plt.plot(years,amipi[e,:],linestyle='-',color='w',linewidth=0.2,
                     alpha=0.6,clip_on=False)
            corrs[e] = sts.pearsonr(indexo,amipi[e])[0]
        print(np.argmax(corrs),'<------ ensemble with highest correlation index!')
            
        plt.plot(years,np.nanmean(amipi,axis=0),linestyle='--',dashes=(1,0.3),color='crimson',
                 linewidth=3.3,clip_on=False,
                 label=r'\textbf{SPEAR_MED_AMIP}')
        plt.plot(years,indexo,linestyle='-',color='deepskyblue',linewidth=1.7,clip_on=False,
                 marker='o',label=r'\textbf{ERA5}',markersize=8)
        
        l = plt.legend(shadow=False,fontsize=12,loc='upper center',
                    fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.515,0.12),
                    labelspacing=1,columnspacing=1,handletextpad=0.5)
        for text in l.get_texts():
            text.set_color("w")
        
        plt.xticks(np.arange(1980,2021,5),map(str,np.arange(1980,2021,5)),fontsize=7,color='w')
        plt.yticks(np.arange(-5,5.5,1),map(str,np.arange(-5,5.5,1)),fontsize=7,color='k')
        plt.xlim([1979,2020])
        plt.ylim([-4,4])
        
        meancorr = np.round(np.mean(corrs),2)
        maxcorr = np.round(np.max(corrs),2)
        mincorr = np.round(np.min(corrs),2)
        plt.text(1979,-4,r'R=%s [%s,%s]' % (meancorr,mincorr,maxcorr),fontsize=8,color='k')
        
        plt.ylabel(r'\textbf{Standardized Index -- T2M}',fontsize=8,color='k')
        
        plt.tight_layout()
        plt.savefig(directoryfigure + 'POSTER_Index-%s_Model-%s_%s_trendIncluded.png' % (indicesdir[i],models[m],monq),dpi=600)


