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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/CoupledSPEAR/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/'

### Loop through indices
indices = ['NINO34']
indicesvari = ['SST']
modelname = ['SPEAR MED']

### Parameters
yearmin = 1979
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)
mon_index = 'JFM'
mon_indexComposite = 'April'
mon_t2m = 'April'
iii = 0 

### Read in SPEAR
amip = np.load(directorydata2 + 'coupledSPEAR/ClimateIndices-Tropical_coupledSPEAR_%s_1979-2019.npz' % mon_index)
indices = amip['indices']
indicesvari = amip['indicesvari']
yearsamip = amip['years']
spear = amip['spear'][iii]

amipComposite  = np.load(directorydata2 + 'coupledSPEAR/ClimateIndices_coupledSPEAR_%s_1979-2019.npz' % mon_indexComposite)
indicesComposite = amipComposite['indices']
indicesvariComposite = amipComposite['indicesvari']
yearsamipComposite = amipComposite['years']
spearComposite = amipComposite['spear'][0]

amiptas = np.load(directorydata2 + 'T2M_BoxNA_AllModels/ClimateIndices_coupledSPEAR-T2M_BoxNA_%s_1979-2019.npz' % mon_index)
speartas = amiptas['spear'][0]

###############################################################################
###############################################################################
###############################################################################
### Read in SPEAR-NOAER
amipaer = np.load(directorydata2 + 'coupledSPEAR-NOAER/ClimateIndices-Tropical_coupledSPEAR-NOAER_%s_1979-2019.npz' % mon_t2m)
indicesaer = amipaer['indices']
indicesvariaer = amipaer['indicesvari']
yearsamipaer = amipaer['years']
spearaer = amipaer['spear'][iii]

amipaerComposite  = np.load(directorydata2 + 'coupledSPEAR-NOAER/ClimateIndices_coupledSPEAR-NOAER_%s_1979-2019.npz' % mon_indexComposite)
indicesaerComposite = amipaerComposite['indices']
indicesaervariComposite = amipaerComposite['indicesvari']
yearsaeramipComposite = amipaerComposite['years']
spearaerComposite = amipaerComposite['spear'][0]

amipaertas = np.load(directorydata2 + 'T2M_BoxNA_AllModels/ClimateIndices_coupledSPEAR-NOAER-T2M_BoxNA_%s_1979-2019.npz' % mon_index)
spearaertas = amipaertas['spear'][0]

def compositePhase(data,index,indexComposite,splitmin,splitmax):
    posneg = []
    negpos = []
    for e in range(data.shape[0]):
        for yr in range(data.shape[1]):
            if index[e,yr] > splitmax:
                if indexComposite[e,yr] < splitmax:
                    posneg.append(data[e,yr])
            elif index[e,yr] < splitmin:
                if indexComposite[e,yr] > splitmax:
                    negpos.append(data[e,yr])
    return posneg,negpos

### Composite onto indices
splitmin = -0.5
splitmax = 0.5
spear_posneg,spear_negpos = compositePhase(speartas,spear,spearComposite,splitmin,splitmax)
spearaer_posneg,spearaer_negpos = compositePhase(spearaertas,spearaer,spearaerComposite,splitmin,splitmax)

def calcPDFs(pos,neg):
    num_bins = np.arange(-4,4.01,0.01)
    
    m_pos,s_pos = sts.norm.fit(pos)
    pdf_pos = sts.norm.pdf(num_bins,m_pos,s_pos)
    
    m_neg,s_neg = sts.norm.fit(neg)
    pdf_neg = sts.norm.pdf(num_bins,m_neg,s_neg)
    
    return pdf_pos,pdf_neg,num_bins
                
### Calculate PDFs
spear_pdf_pos,spear_pdf_neg,bins  = calcPDFs(spear_posneg,spear_negpos)
spearaer_pdf_pos,spearaer_pdf_neg,bins  = calcPDFs(spearaer_posneg,spearaer_negpos)

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
  
plt.plot(bins,spear_pdf_pos,color='crimson',linewidth=2.1,
         linestyle='-',label=r'\textbf{SPEAR-MED [+%s]}' % indices[iii],clip_on=False)
plt.plot(bins,spear_pdf_neg,color='crimson',linewidth=1.4,
         linestyle='--',dashes=(1,0.5),label=r'\textbf{SPEAR-MED [-%s]}' % indices[iii],clip_on=False)  
  
plt.plot(bins,spearaer_pdf_pos,color='deepskyblue',linewidth=2.1,
         linestyle='-',label=r'\textbf{SPEAR-MED-NOAER [+%s]}' % indices[iii],clip_on=False)
plt.plot(bins,spearaer_pdf_neg,color='deepskyblue',linewidth=1.4,
         linestyle='--',dashes=(1,0.5),label=r'\textbf{SPEAR-MED-NOAER [-%s]}' % indices[iii],clip_on=False)

plt.yticks(np.arange(0,1,0.1),list(map(str,np.round(np.arange(0,1,0.1),2))),
           fontsize=6)
plt.xticks(np.arange(-5,5.1,0.5),list(map(str,np.arange(-5,5.1,0.5))),
           fontsize=6) 
plt.xlim([-4,4])
plt.ylim([0,0.5])

l = plt.legend(shadow=False,fontsize=10,loc='upper center',
           fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0.5,1.17),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
for text in l.get_texts():
    text.set_color('k')

plt.ylabel(r'\textbf{Density}',color='k',fontsize=8)  
plt.xlabel(r'\textbf{T2M-BoxNA [$^{\circ}$C]}',color='k',fontsize=8)

plt.savefig(directoryfigure + 'PDFs_%s-Index-%s_T2M-%s_AL-%s_SPEARMED_Combined.png' % (indices[iii],mon_index,mon_t2m,mon_indexComposite),dpi=300)
