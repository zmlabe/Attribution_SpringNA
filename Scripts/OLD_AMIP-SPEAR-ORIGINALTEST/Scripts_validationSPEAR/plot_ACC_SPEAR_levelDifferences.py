"""
Calculate ACC for SPEAR_33lev
 
Author    : Zachary M. Labe
Date      : 30 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import calc_DetrendData as DT
import sys
import itertools
import read_ERA5_monthlyLOWS_hindcastSPEAR_65lev as ER
import read_SPEAR_LOW_Hindcast_33lev as SPlow
import read_SPEAR_LOW_Hindcast_65lev as SPhigh
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/SPEAR_Low_65lev/' 
directorydata = '/work/Zachary.Labe/Data/SPEAR/Hindcasts/SPEAR_LOW/' 

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

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
experi1 = ['SPEAR-LOW-65Lev','SPEAR-LOW-33Lev','Difference']
variq = 'T2M'
slicemonthall = ['FMA','MA','April']
# slicemonthall = ['MA','April']
monthhind = 'Jan'
slicenan = 'nan'
detrend = False
yearmin = 1995
yearmax = 2019
years = np.arange(1995,2019+1,1)

### Read in obs
for yr in range(len(slicemonthall)):
    slicemonth = slicemonthall[yr]
    latobs,lonobs,levobs,obsq = ER.read_ERA5LOWS_monthlyHindcast(variq,'/work/Zachary.Labe/Data/',monthhind,slicenan,'surface')
    
    ### Calculate anomalies
    clim = np.nanmean(obsq,axis=0)
    obs = obsq - clim
        
    ### Calculate time periods
    if monthhind == 'Jan':
        if slicemonth == 'FMA':
            obs_time = np.nanmean(obs[:,1:4,:,:],axis=1)
        elif slicemonth == 'MA':
            obs_time = np.nanmean(obs[:,2:4,:,:],axis=1)
        elif slicemonth == 'April':
            obs_time = np.nanmean(obs[:,3:4,:,:],axis=1)
    elif monthhind == 'Feb':
        if slicemonth == 'FMA':
            obs_time = np.nanmean(obs[:,0:3,:,:],axis=1)
        elif slicemonth == 'MA':
            obs_time = np.nanmean(obs[:,1:3,:,:],axis=1)
        elif slicemonth == 'April':
            obs_time = np.nanmean(obs[:,2:3,:,:],axis=1)
    elif monthhind == 'Mar':
        if slicemonth == 'FMA':
            print(ValueError('FULL TIME PERIOD NOT AVAILABLE!'))
            sys.exit()
        elif slicemonth == 'MA':
            obs_time = np.nanmean(obs[:,0:2,:,:],axis=1)
        elif slicemonth == 'April':
            obs_time = np.nanmean(obs[:,1:2,:,:],axis=1)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################   
    ### Read in initialized run
    lat,lon,spearlow = SPlow.read_SPEAR_LOW_Hindcast33lev(directorydata,variq,monthhind,slicenan,15)
    lat,lon,spearhigh = SPhigh.read_SPEAR_LOW_Hindcast65lev(directorydata,variq,monthhind,slicenan,15)
    
    def calc_Corrs(spear,lat,lon,obs):
        ### Calculate time periods
        if monthhind == 'Jan':
            if slicemonth == 'FMA':
                spear_time = np.nanmean(spear[:,:,1:4,:,:],axis=2)
            elif slicemonth == 'MA':
                spear_time = np.nanmean(spear[:,:,2:4,:,:],axis=2)
            elif slicemonth == 'April':
                spear_time = np.nanmean(spear[:,:,3:4,:,:],axis=2)
        elif monthhind == 'Feb':
            if slicemonth == 'FMA':
                spear_time = np.nanmean(spear[:,:,0:3,:,:],axis=2)
            elif slicemonth == 'MA':
                spear_time = np.nanmean(spear[:,:,1:3,:,:],axis=2)
            elif slicemonth == 'April':
                spear_time = np.nanmean(spear[:,:,2:3,:,:],axis=2)
        elif monthhind == 'Mar':
            if slicemonth == 'FMA':
                print(ValueError('FULL TIME PERIOD NOT AVAILABLE!'))
                sys.exit()
            elif slicemonth == 'MA':
                spear_time = np.nanmean(spear[:,:,0:2,:,:],axis=2)
            elif slicemonth == 'April':
                spear_time = np.nanmean(spear[:,:,1:2,:,:],axis=2)
                
        ### Calculate means across the ensembles
        spear_hind = np.nanmean(spear_time[:,:,:,:],axis=1)
        
        ### Regrid obs for comparison
        lon2,lat2 = np.meshgrid(lon,lat)
        newobs = obs_time
        
        if detrend == True:
            obsdt = DT.detrendDataR(newobs,'surface','monthly')
            moddt = DT.detrendDataR(spear_hind,'surface','monthly')
        else:
            obsdt = newobs
            moddt = spear_hind
            
        ### Begin function to correlate observations with model, ensemble, year
        corr = np.empty((moddt.shape[1],moddt.shape[2]))
        pvals = np.empty((moddt.shape[1],moddt.shape[2]))
        for i in range(moddt.shape[1]):
            for j in range(moddt.shape[2]):
                xx = obsdt[:,i,j]
                yy = moddt[:,i,j]
                if np.isfinite(np.nanmax(xx)) == True and np.isfinite(np.nanmax(yy)) == True:
                    na = np.logical_or(np.isnan(xx),np.isnan(yy))
                    corr[i,j],pvals[i,j]=sts.pearsonr(xx[~na],yy[~na])
                else:
                    corr[i,j] = np.nan
                    pvals[i,j] = np.nan
                
        return corr,pvals,lat2,lon2
    
    ### Return ACC
    corr_low,pvals_low,lat2,lon2 = calc_Corrs(spearlow,lat,lon,obs_time)
    corr_high,pvals_high,lat2,lon2 = calc_Corrs(spearhigh,lat,lon,obs_time)
    diff = corr_high - corr_low
            
    varplot = [corr_high,corr_low,diff]
    varp = [pvals_high,pvals_low]

    #######################################################################
    #######################################################################
    #######################################################################
    ### Plot subplot of mean climate
    if variq == 'T2M':
        limit = np.arange(-1,1.01,0.1)
        barlim = np.round(np.arange(-1,1.1,0.5),2)
        cmap = cmocean.cm.balance
        label = r'\textbf{Anomaly Correlation Coefficient}'
        
    fig = plt.figure(figsize=(8,4))
    for pl in range(len(varplot)):
        var = varplot[pl]
    
        ax1 = plt.subplot(1,3,pl+1)
       
        m = Basemap(projection='ortho',lon_0=250,lat_0=4z0,resolution='l',area_thresh=10000)
        m.drawcoastlines(color='dimgrey',linewidth=1)
        m.drawstates(color='dimgrey',linewidth=0.5)
        m.drawcountries(color='dimgrey',linewidth=0.5)
        
        var, lons_cyclic = addcyclic(var, lon)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
                  
        circle = m.drawmapboundary(fill_color='white',color='dimgrey',
                                   linewidth=0.7)
        circle.set_clip_on(False)
        
        cs1 = m.contourf(x,y,var,limit,extend='both')
        
        if pl < 2:
            pvar = varp[pl]
            ### Significant at 95% confidence level
            pvar[np.where(pvar >= 0.05)] = np.nan
            pvar[np.where(pvar < 0.05)] = 1.
            pvar,lons_cyclic = addcyclic(pvar, lon)
            csp = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                         linewidths=0.4)
        
        cs1.set_cmap(cmocean.cm.balance)
        
        ## Box 1
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
        latsslice = np.ones(len(lonsslice))*la2
        m.plot(lonsslice, latsslice, color='gold', linewidth=2, latlon=True,zorder=4)
        latsslice = np.ones(len(lonsslice))*la1
        m.plot(lonsslice, latsslice, color='gold', linewidth=2, latlon=True,zorder=4)
        m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=2,color='gold',zorder=4)
        m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=2,color='gold',zorder=4)
        
        plt.title(r'%s' % experi1[pl],fontsize=12,color='dimgrey')
        ax1.annotate(r'\textbf{[%s]}' % (letters[pl]),xy=(0,0),xytext=(0.07,0.90),
                  textcoords='axes fraction',color='k',fontsize=9,
                  rotation=40,ha='center',va='center')
            
    fig.suptitle(r'\textbf{%s; initialized on %s}' % (slicemonth,monthhind),fontsize=18,color='k')
    
    ###############################################################################
    cbar_ax1 = fig.add_axes([0.41,0.09,0.2,0.035])                
    cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                        extend='neither',extendfrac=0.07,drawedges=False)
    cbar1.set_label(label,fontsize=9,color='k',labelpad=1.4)  
    cbar1.set_ticks(barlim)
    cbar1.set_ticklabels(list(map(str,barlim)))
    cbar1.ax.tick_params(axis='x', size=.001,labelsize=7,pad=2)
    cbar1.outline.set_edgecolor('dimgrey')
    
    plt.tight_layout()
    if detrend == False:
        plt.savefig(directoryfigure + 'ACC_%s_levelDifferences_init-%s_%s_%s-%s.png'% (slicemonth,monthhind,variq,yearmin,yearmax),dpi=300)
    else: 
        plt.savefig(directoryfigure + 'ACC_%s_levelDifferences_init-%s_%s_%s-%s_detrend.png'% (slicemonth,monthhind,variq,yearmin,yearmax),dpi=300)
