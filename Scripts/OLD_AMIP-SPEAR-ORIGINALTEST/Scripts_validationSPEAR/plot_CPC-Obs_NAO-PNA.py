"""
Plot monthly trends in the Arctic

Author    : Zachary M. Labe
Date      : 27 June 022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/' 

monq = 'FMA'
yearsall = np.arange(1950,2021+1,1)
yearmin = 1979
yearmax = 2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

### Read in observations
obs = np.genfromtxt(directorydata2 + 'ERA5_T2M_BoxOfInterest_TimeSeries-%s.txt' % monq,unpack=True)
obsz = sts.zscore(obs)

### Read in teleconnections
aod = np.genfromtxt(directorydata + 'AO/AO_CPC_1950-2021.txt',unpack=True)
years_ao = aod[0,:]
ao_index = aod[1:,yearq]

naod = np.genfromtxt(directorydata + 'NAO/NAO_CPC_1950-2021.txt',unpack=True)
years_nao = naod[0,:]
nao_index = naod[1:,yearq]

pnad = np.genfromtxt(directorydata + 'PNA/PNA_CPC_1950-2021.txt',unpack=True)
years_pna = pnad[0,:]
pna_index = pnad[1:,yearq]

### Slice months
if monq == 'FMA':
    ao = np.nanmean(ao_index[1:4,:],axis=0)
    nao = np.nanmean(nao_index[1:4,:],axis=0)
    pna = np.nanmean(pna_index[1:4,:],axis=0)
elif monq == 'April':
    ao = np.nanmean(ao_index[3:4,:],axis=0)
    pna = np.nanmean(pna_index[3:4,:],axis=0)
    nao = np.nanmean(nao_index[3:4,:],axis=0)
elif monq == 'none':
    ao = ao_index
    pna = pna_index
    nao = nao_index
else:
    print(ValueError('wrong months selected!'))
    sys.exit()
    
### Calculate trends
slope_ao,intercept_ao,r_ao,p_ao,se_ao = sts.linregress(years,ao)
slope_nao,intercept_nao,r_nao,p_nao,se_nao = sts.linregress(years,nao)
slope_pna,intercept_pna,r_pna,p_pna,se_pna = sts.linregress(years,pna)

### Calculate correlations
corr_ao,pcorr_ao = sts.pearsonr(obsz,ao)
corr_nao,pcorr_nao = sts.pearsonr(obsz,nao)
corr_pna,pcorr_pna = sts.pearsonr(obsz,pna)
    
### Calculate class combinations
classes_AO = np.empty((years.shape[0]))
for i in range(years.shape[0]):
    if ((ao[i] >= 0.05) and (pna[i] >= 0.05)):
        classes_AO[i] = 0
    elif ((ao[i] <= -0.05) and (pna[i] <= -0.05)):
        classes_AO[i] = 1
    elif ((ao[i] >= -0.05) and (pna[i] <= -0.05)):
        classes_AO[i] = 2
    elif ((ao[i] <= -0.05) and (pna[i] >= -0.05)):
        classes_AO[i] = 3
    else:
        classes_AO[i] = 4
        
classes_NAO = np.empty((years.shape[0]))
for i in range(years.shape[0]):
    if ((nao[i] >= 0.05) and (pna[i] >= 0.05)):
        classes_NAO[i] = 0
    elif ((nao[i] <= -0.05) and (pna[i] <= -0.05)):
        classes_NAO[i] = 1
    elif ((nao[i] >= -0.05) and (pna[i] <= -0.05)):
        classes_NAO[i] = 2
    elif ((nao[i] <= -0.05) and (pna[i] >= -0.05)):
        classes_NAO[i] = 3
    else:
        classes_NAO[i] = 4
        
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

rects = plt.bar(years,pna,align='center',zorder=2,color='crimson',label=r'PNA',clip_on=False)
plt.plot(years,ao,linestyle='-',color='darkblue',linewidth=1,clip_on=False,label='AO')
plt.plot(years,nao,linestyle='--',dashes=(1,0.3),color='deepskyblue',linewidth=1,clip_on=False,label='NAO')
plt.plot(years,obsz,linestyle='--',dashes=(1,0.3),color='k',linewidth=2,clip_on=False,label='ERA5-T2M')

l = plt.legend(shadow=False,fontsize=7,loc='upper center',
            fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,0.07),
            labelspacing=1,columnspacing=1,handletextpad=0.4)

plt.xticks(np.arange(1980,2021,5),map(str,np.arange(1980,2021,5)),fontsize=7)
plt.yticks(np.arange(-5,5.5,0.5),map(str,np.arange(-5,5.5,0.5)),fontsize=7)
plt.xlim([1979,2020])
plt.ylim([-3,3])

plt.ylabel(r'\textbf{Standardized Index}',fontsize=8,color='k')
plt.title(r'\textbf{Observational Indices For %s}' % monq,fontsize=12,color='k')
plt.savefig(directoryfigure + 'CPC_Indices-Obs_%s' % monq,dpi=300)


