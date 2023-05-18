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
import cmasher as cmr

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/Presentations/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 

monq = 'April'
yearsall = np.arange(1950,2021+1,1)
yearmin = 1979
yearmax = 2020
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

### Read in observations
obs = np.genfromtxt(directorydata + 'T2M_BoxNA/T2M_BoxNA_%s_%s_%s-%s_detrended.txt' % ('T2M',monq,
                                                                yearmin,yearmax),
                                                                unpack=True)
obsz = sts.zscore(obs[1,:])

### Read in indices
indexallname = ['T2M-Reg','NECPSST','WPPRECT','NINO34','NPI','PCHT2M','PCHZ100','PCHZ30','PCHZ1000','SHI','UBI','AL','PNA','NAO','NPO','ABI','NPOcustomZ']
variqn = ['SST','P','SST','SLP','T2M','Z100','Z30','Z1000','SLP','Z500','SLP','Z500','SLP','SLP','Z500','Z500']

indexall = np.empty((len(indexallname)-1,years.shape[0]))
for i in range(len(indexallname)-1):
    if any([indexallname[i+1]=='PNA',indexallname[i+1]=='NAO']):
        indexq = np.genfromtxt(directorydata + '%s/%smodified_%s_%s_%s-%s_detrended.txt' % (indexallname[i+1],
                                                                        indexallname[i+1],variqn[i],monq,
                                                                        yearmin,yearmax),
                                                                        unpack=True)
    else:
        indexq = np.genfromtxt(directorydata + '%s/%s_%s_%s_%s-%s_detrended.txt' % (indexallname[i+1],
                                                                        indexallname[i+1],variqn[i],monq,
                                                                        yearmin,yearmax),
                                                                        unpack=True)
    indexall[i,:] = sts.zscore(indexq[1,:])
    print('Read in data for ---> %s' % indexallname[i+1])
    
### All climate indices
tempindexcpc = np.concatenate([obsz[np.newaxis,:],indexall],axis=0)

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(111)

indexallnamen = ['T2M-Reg','NECP SST','WP PRECT','NINO34','NPI','T2M-Arctic','PCHZ100','PCHZ30','PCHZ1000','SHI','UBI','AL','PNA','NAO','NPO','ABI','NPO*']

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

cs = plt.pcolormesh(tempindexcpc,shading='faceted',edgecolor='k',
                    linewidth=0.3,vmin=-3,vmax=3)

# for i in range(tempindexcpc.shape[0]):
#     for j in range(tempindexcpc.shape[1]):
#         plt.text(j+0.5,i+0.5,r'\textbf{%+1.2f}' % tempindexcpc[i,j],fontsize=7,
#                   color='k',va='center',ha='center')

cs.set_cmap(cmr.fusion_r)

ylabels = indexallnamen
plt.yticks(np.arange(0.5,len(indexallnamen)+0.5,1),ylabels,ha='right',color='w',
            va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)
plt.xticks(np.arange(0.5,42.5,4),np.arange(1979,2019+1,4),ha='center',color='w',
            va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,41])

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
ticks = np.arange(-3,4,1)
labels = list(map(str,np.arange(-3,4,1)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('darkgrey')
cbar.set_label(r'\textbf{Standardized Climate Indices - %s}' % monq,
                color='darkgrey',labelpad=3,fontsize=12)

plt.tight_layout()
plt.savefig(directoryfigure + 'PRESENT_HeatMap_ClimateIndices-Obs_%s' % monq,dpi=900)

###############################################################################
###############################################################################
###############################################################################
### Calculate and plot correlation matrix
yearmn = 1979
yearmx = 2020
yearq = np.where((years >= yearmn) & (years <= yearmx))[0]
# yearq = np.where(years==2020)[0][0] + 1
corr=np.corrcoef(tempindexcpc[:,yearq])
mask = np.zeros_like(corr,dtype=np.bool).T
corr[np.triu_indices_from(mask)] = np.nan

fig = plt.figure(figsize=(8,4))
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

cs = plt.pcolormesh(corr,shading='faceted',edgecolor='k',
                    linewidth=1.7,vmin=-1,vmax=1)

for i in range(corr.shape[0]):
    for j in range(corr.shape[1]):
        if np.isnan(corr[i,j]) == False:
            if (corr[i,j] > 0.305) or (corr[i,j] <-0.305):
                print(corr[i,j])
                if corr[i,j] > 0:
                    cc = 'k'
                else:
                    cc = 'k'
                plt.text(j+0.5,i+0.5-0.1,r'\textbf{%+1.2f}' % corr[i,j],fontsize=8,
                         color=cc,va='center',ha='center')

cs.set_cmap(cmr.fusion_r)

ylabels = indexallnamen
xlabels = indexallnamen
plt.yticks(np.arange(0.5,len(indexallnamen)+0.5,1),ylabels,ha='right',color='w',
           va='center',size=7)
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)

plt.xticks(np.arange(0.5,len(indexallnamen)+0.5,1),xlabels,ha='center',color='w',
           va='top',rotation=90,size=7)
xax = ax.get_xaxis()
xax.set_tick_params(pad=2)

cbar = plt.colorbar(cs,orientation='horizontal',shrink = 0.5)
ticks = np.arange(-1,1.1,0.2)
labels = list(map(str,np.round(np.arange(-1,1.1,0.2),2)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels,fontsize=7)
cbar.ax.tick_params(axis='x',size=.001)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{Correlation Matrix -- April}',
               color='darkgrey',labelpad=3,fontsize=12)

plt.tight_layout()
plt.savefig(directoryfigure + 'PRESENT_CorrelationMatrix_ClimateIndices-Obs_%s_lead_%s-%s' % (monq,yearmn,yearmx),dpi=300)
