"""
Plot heatmap of climate indices for 1979-2019

Author    : Zachary M. Labe
Date      : 1 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts
import cmocean

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 

monq = 'April'
yearsall = np.arange(1950,2021+1,1)
yearmin = 1979
yearmax = 2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

### Read in observations
obs = np.genfromtxt(directorydata + 'T2M_BoxNA/ERA5_T2M_BoxOfInterest_TimeSeries-%s.txt' % 'FMA',unpack=True)
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

### Read in indices
indexallname = ['T2M-Reg','NECPSST','WPPRECT','NINO34','NPI','PCHT2M','PCHZ100','PCHZ30','PCHZ1000','SHI','UBI','AL','PNA','NAO','AO']
variqn = ['SST','P','SST','SLP','T2M','Z100','Z30','Z1000','SLP','Z500','SLP']

indexall = np.empty((len(indexallname)-4,years.shape[0]))
for i in range(len(indexallname)-4):
    indexq = np.genfromtxt(directorydata + '%s/%s_%s_%s_%s-%s_detrended.txt' % (indexallname[i+1],
                                                                    indexallname[i+1],variqn[i],monq,
                                                                    yearmin,yearmax),
                                                                    unpack=True)
    indexall[i,:] = sts.zscore(indexq[1,:])
    print('Read in data for ---> %s' % indexallname[i+1])
    
### All climate indices
tempindexcpc = np.concatenate([obsz[np.newaxis,:],indexall,pna[np.newaxis,:],nao[np.newaxis,:],ao[np.newaxis,:]],axis=0)

# ###############################################################################
# ###############################################################################
# ###############################################################################
# fig = plt.figure()
# ax = plt.subplot(111)

# ax.spines['top'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
#             color='w')

# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=True,      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom=True)
# plt.tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     left=False,      # ticks along the bottom edge are off
#     right=False,         # ticks along the top edge are off
#     labelleft='on')

# cs = plt.pcolormesh(tempindexcpc,shading='faceted',edgecolor='w',
#                     linewidth=0.3,vmin=-3,vmax=3)

# # for i in range(tempindexcpc.shape[0]):
# #     for j in range(tempindexcpc.shape[1]):
# #         plt.text(j+0.5,i+0.5,r'\textbf{%+1.2f}' % tempindexcpc[i,j],fontsize=7,
# #                  color='k',va='center',ha='center')

# cs.set_cmap(cmocean.cm.balance)

# ylabels = indexallname
# plt.yticks(np.arange(0.5,len(indexallname)+0.5,1),ylabels,ha='right',color='dimgrey',
#            va='center')
# yax = ax.get_yaxis()
# yax.set_tick_params(pad=-2)
# plt.xticks(np.arange(0.5,42.5,4),np.arange(1979,2019+1,4),ha='center',color='dimgrey',
#            va='center')
# xax = ax.get_xaxis()
# xax.set_tick_params(pad=8)
# plt.xlim([0,41])

# cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
# ticks = np.arange(-3,4,1)
# labels = list(map(str,np.arange(-3,4,1)))
# cbar.set_ticks(ticks)
# cbar.set_ticklabels(labels)
# cbar.ax.tick_params(axis='x', size=.001)
# cbar.outline.set_edgecolor('dimgrey')
# cbar.set_label(r'\textbf{Standardized Climate Indices - %s}' % monq,
#                color='dimgrey',labelpad=3,fontsize=12)

# plt.tight_layout()
# plt.savefig(directoryfigure + 'HeatMap_ClimateIndices-Obs_%s' % monq,dpi=900)

###############################################################################
###############################################################################
###############################################################################
### Calculate and plot correlation matrix
yearq = np.where(years==2013)[0][0] + 1
corr=np.corrcoef(tempindexcpc[:,:yearq])
mask = np.zeros_like(corr,dtype=np.bool).T
corr[np.triu_indices_from(mask)] = np.nan

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

cs = plt.pcolormesh(corr,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-1,vmax=1)

for i in range(corr.shape[0]):
    for j in range(corr.shape[1]):
        if np.isnan(corr[i,j]) == False:
            plt.text(j+0.5,i+0.5,r'\textbf{%+1.2f}' % corr[i,j],fontsize=6,
                     color='dimgrey',va='center',ha='center')

cs.set_cmap(cmocean.cm.balance)

ylabels = indexallname
xlabels = indexallname
plt.yticks(np.arange(0.5,len(indexallname)+0.5,1),ylabels,ha='right',color='dimgrey',
           va='center',size=7)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)

plt.xticks(np.arange(0.5,len(indexallname)+0.5,1),xlabels,ha='center',color='dimgrey',
           va='top',rotation=90,size=7)
xax = ax.get_xaxis()
xax.set_tick_params(pad=2)

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
ticks = np.arange(-1,1.1,0.2)
labels = list(map(str,np.round(np.arange(-1,1.1,0.2),2)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels,fontsize=7)
cbar.ax.tick_params(axis='x',size=.001)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{Correlation Matrix - %s for April-T2M}' % monq,
               color='dimgrey',labelpad=3,fontsize=12)

plt.tight_layout()
plt.savefig(directoryfigure + 'CorrelationMatrix_ClimateIndices-Obs_%s_lead_%s' % (monq,yearq-1),dpi=900)
