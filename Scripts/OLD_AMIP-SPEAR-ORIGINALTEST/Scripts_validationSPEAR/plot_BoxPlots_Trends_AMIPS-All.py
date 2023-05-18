"""
Plot box plots of trends over different time periods compared to SMILEs

Author    : Zachary M. Labe
Date      : 23 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import scipy.stats as sts
import itertools

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/'

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
scenario = ['amip_1880s_rf','amip_obs_rf','amip_clim_polar','spear']
model_1880s_rf = ['CAM4','ECHAM5','ESRL-CAM5']
model_obs_rf = ['CAM4','ECHAM5','ESRL-CAM5']
model_clim_polar = ['CAM4','ECHAM5','ESRL-CAM5','ESRL-GFS']
model_spear = ['spear']
variq = 'T2M'
slicemonth = 'FMA'
model = [model_1880s_rf,model_obs_rf,model_clim_polar,model_spear]
modelunravel = list(itertools.chain(*model))
yearsobs_all = np.arange(1979,2014+1,1)
yearsobs_rec = np.arange(2000,2014+1,1)

### Read in data
trend_obs_all = np.loadtxt(directoryoutput + 'Slopes_%s_obs_%s-%s.txt' % (slicemonth,1979,2014))
trend_obs_rec = np.loadtxt(directoryoutput + 'Slopes_%s_obs_%s-%s.txt' % (slicemonth,2000,2014))
amipdata = np.load(directoryoutput + 'AMIPs/AMIP_Trends_%s_%s.npz' % (slicemonth,variq),
                   allow_pickle=True)
year = np.load(directoryoutput + 'AMIPs/AMIP_Years_%s_%s.npy' % (slicemonth,variq),
                   allow_pickle=True)

trend_all = amipdata['trend'][0]
intercept_all = amipdata['intercept'][0]
trend_rec = amipdata['trend'][1]
intercept_rec = amipdata['intercept'][1]

decobs_all = trend_obs_all*10.
decobs_rec = trend_obs_rec*10.

### Plot max year
model_1880s_rfy = ['CAM4(%s)'% (year[0][0].max()),
                   'ECHAM5(%s)'% (year[0][1].max()),
                   'CAM5(%s)'% (year[0][2].max())]
model_obs_rfy = ['CAM4(%s)'% (year[1][0].max()),
                 'ECHAM5(%s)'% (year[1][1].max()),
                 'CAM5(%s)'% (year[1][2].max())]
model_clim_polary = ['CAM4(%s)'% (year[2][0].max()),
                     'ECHAM5(%s)'% (year[2][1].max()),
                     'CAM5(%s)'% (year[2][2].max()),
                     'GFS(%s)'% (year[2][3].max())]
model_speary = ['spear(%s)'% (year[3][0].max())]
modely = [model_1880s_rfy,model_obs_rfy,model_clim_polary,model_speary]

###############################################################################
###############################################################################
###############################################################################               
### Plot Figure
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

fig = plt.figure(figsize=((10,5)))
for vv in range(len(scenario)-1):
    if vv < 2:
        axb = plt.subplot(1,3,vv+1)
        
        axb.spines['top'].set_color('none')
        axb.spines['right'].set_color('none')
        axb.spines['left'].set_color('none')
        axb.spines['bottom'].set_color('none')
        axb.tick_params('both',length=4,width=1,which='major',color='darkgrey')
        axb.xaxis.grid(zorder=1,color='darkgrey',alpha=1,linewidth=1,clip_on=False)
        plt.axhline(y=0,color='darkgrey',linestyle='-',linewidth=1,zorder=1)
        
        datum = [decade*10 for decade in trend_all[vv]]
        
        vp = plt.violinplot(datum,showmeans=False,showmedians=True,
                            vert=True,widths=0.6,showextrema=True)
        
        if vv == 1:
            spear = [decade*10 for decade in trend_all[-1]]
            for spr in range(len(spear[0])):
                plt.axhline(spear[0][spr],color='b',linestyle='-',
                            linewidth=1)
        
        plt.axhline(decobs_all,color='k',linestyle='--',dashes=(1,0.3),linewidth=2)
        plt.setp(axb,xticks=[x+1 for x in range(len(datum))],
                              xticklabels=modely[vv])
        plt.setp(axb.get_xticklabels(),fontsize=7)
        
        positionsq = np.array(np.arange(1,4,1))
        for i in range(len(datum)):
            x = datum[i]
            y = np.random.normal(positionsq[i], 0.02, size=len(x))
            plt.plot(y,x,color='k',alpha=0.5,zorder=10,marker='.',linewidth=0,markersize=10,markeredgewidth=0,clip_on=False)
                             
        for i in vp['bodies']:
            i.set_edgecolor('dimgrey')  
        vp['cbars'].set_color('k')
        vp['cmaxes'].set_color('k')
        vp['cmins'].set_color('k')
        vp['cmedians'].set_color('k')
        vp['cmedians'].set_linewidth(3)
        vp['cmaxes'].set_linewidth(3)        
        vp['cmins'].set_linewidth(3)       
        vp['cmaxes'].set_linestyle('-')        
        vp['cmins'].set_linestyle('-')          
        vp['bodies'][0].set_facecolor('maroon')
        vp['bodies'][1].set_facecolor('teal')   
        vp['bodies'][2].set_facecolor('goldenrod')      
        vp['bodies'][0].set_alpha(0.8)
        vp['bodies'][1].set_alpha(0.8)  
        vp['bodies'][2].set_alpha(0.8)  
        
        plt.yticks(np.arange(-2,2.5,0.5),map(str,np.round(np.arange(-2,2.5,0.5),2)),size=7)
        plt.ylim([-1,1])  
        
        if vv == 0:
            plt.ylabel(r'\textbf{%s for %s of T2M Trends [$^{\circ}$C/Decade]' % (slicemonth,scenario[vv]),color='k',size=6)
        
        axb.xaxis.set_ticks_position('bottom')
        axb.yaxis.set_ticks_position('left')
        
        plt.title(r'\textbf{%s}' % scenario[vv],color='k',fontsize=13)
        
    if vv == 2:
        axb = plt.subplot(1,3,vv+1)
        
        axb.spines['top'].set_color('none')
        axb.spines['right'].set_color('none')
        axb.spines['left'].set_color('none')
        axb.spines['bottom'].set_color('none')
        axb.tick_params('both',length=4,width=1,which='major',color='darkgrey')
        axb.xaxis.grid(zorder=1,color='darkgrey',alpha=1,linewidth=1,clip_on=False)
        plt.axhline(y=0,color='darkgrey',linestyle='-',linewidth=1,zorder=1)
        
        datum = [decade*10 for decade in trend_all[vv]]
        
        vp = plt.violinplot(datum,showmeans=False,showmedians=True,
                            vert=True,widths=0.6,showextrema=True)
        plt.axhline(decobs_all,color='k',linestyle='--',dashes=(1,0.3),linewidth=2)
        plt.setp(axb,xticks=[x+1 for x in range(len(datum))],
                              xticklabels=modely[vv])
        plt.setp(axb.get_xticklabels(),fontsize=7)
        
        positionsq = np.array(np.arange(1,5,1))
        for i in range(len(datum)):
            x = datum[i]
            y = np.random.normal(positionsq[i], 0.02, size=len(x))
            plt.plot(y,x,color='k',alpha=0.5,zorder=10,marker='.',linewidth=0,markersize=10,markeredgewidth=0,clip_on=False)
                             
        for i in vp['bodies']:
            i.set_edgecolor('dimgrey')  
        vp['cbars'].set_color('k')
        vp['cmaxes'].set_color('k')
        vp['cmins'].set_color('k')
        vp['cmedians'].set_color('k')
        vp['cmedians'].set_linewidth(3)
        vp['cmaxes'].set_linewidth(3)        
        vp['cmins'].set_linewidth(3)       
        vp['cmaxes'].set_linestyle('-')        
        vp['cmins'].set_linestyle('-')          
        vp['bodies'][0].set_facecolor('maroon')
        vp['bodies'][1].set_facecolor('teal')   
        vp['bodies'][2].set_facecolor('goldenrod')  
        vp['bodies'][3].set_facecolor('deepskyblue')  
        vp['bodies'][0].set_alpha(0.8)
        vp['bodies'][1].set_alpha(0.8)  
        vp['bodies'][2].set_alpha(0.8) 
        vp['bodies'][3].set_alpha(0.8) 
        
        plt.yticks(np.arange(-2,2.5,0.5),map(str,np.round(np.arange(-2,2.5,0.5),2)),size=7)
        plt.ylim([-1,1])  
        
        axb.xaxis.set_ticks_position('bottom')
        axb.yaxis.set_ticks_position('left')
        
        plt.title(r'\textbf{%s}' % scenario[vv],color='k',fontsize=13)

plt.tight_layout()
plt.savefig(directoryfigure + 'AMIPS_TrendViolins_%s_All.png' % (slicemonth),dpi=300)







