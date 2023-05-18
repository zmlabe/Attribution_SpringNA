"""
Calculate ACC for SPEAR-AMIP-old
 
Author    : Zachary M. Labe
Date      : 9 June 2022
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
import read_ERA5_monthly1x1 as ER
import read_FACTS_AMIPS as FA
from scipy.interpolate import griddata as g
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/SPEAR_Low_65lev/' 
directorydata = '/work/Zachary.Labe/Data/SPEAR/Hindcasts/SPEAR_LOW/SPEAR-65lev/' 

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

def regrid(lat1,lon1,lat2,lon2,var,years):
    """
    Interpolated on selected grid. Reads NCEP in as 3d with 
    [year,lat,lon]
    """
    
    varn_re = np.reshape(var,(var.shape[0],(lat1.shape[0]*lon1.shape[1])))   
    
    print('Completed: Start regridding process:')
    varn = np.empty((var.shape[0],lat2.shape[0],lon2.shape[1]))
    for i in range(var.shape[0]):
        z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,:],(lat2,lon2),method='linear')
        varn[i,:,:] = z
        print('Completed: Year %s Regridding---' % (years[i]))
    return varn

### Parameters
variq = 'T2M'
slicemonthall = ['FMA','MA','April']
slicenan = 'nan'
detrend = False
yearmin = 1995
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)
yearsoriginal = np.arange(1979,2021+1,1)
yearsorigina2 = np.arange(1921,2020+1,1)

### Read in obs
for yr in range(len(slicemonthall)):
    slicemonth = slicemonthall[yr]
    latobs,lonobs,levobs,var = ER.read_ERA5_monthly1x1(variq,'/work/Zachary.Labe/Data/',slicemonth,
                                    years,3,True,slicenan,'surface')
    
    ### Slice only years from 1979 - 2019
    yearq = np.where((yearsoriginal >= yearmin) & (yearsoriginal <= yearmax))[0]
    obsq = var[yearq,:,:]
    
    ### Calculate anomalies
    clim = np.nanmean(obsq,axis=0)
    obs = obsq - clim
    
    ###############################################################################
    ###############################################################################
    ###############################################################################   
    ### Read in AMIP 
    lat,lon,lev,enseall,yearsf = FA.read_FACTS_Experi('spear','spear',variq,slicemonth,4,slicenan,'surface')
    
    ### Slice only years from 1979 - 2019
    yearq2 = np.where((yearsf >= yearmin) & (yearsf <= yearmax))[0]
    spearpre = enseall[:,yearq2,:,:]
    
    ### Calculate anomalies
    clim = np.nanmean(spearpre,axis=1)[:,np.newaxis,:,:]
    spear = spearpre - clim
    
    ### Regrid obs for comparison
    lon2,lat2 = np.meshgrid(lon,lat)
    lonobs2,latobs2 = np.meshgrid(lonobs,latobs)
    newobs = regrid(latobs2,lonobs2,lat2,lon2,obs,years)  
    
    if detrend == True:
        obsdt = DT.detrendDataR(newobs,'surface','monthly')
        moddt = DT.detrendData(spear,'surface','monthly')
    else:
        obsdt = newobs
        moddt = spear
       
    ### Begin function to correlate observations with model, ensemble, year
    corr = np.empty((moddt.shape[0],moddt.shape[2],moddt.shape[3]))
    pvals = np.empty((moddt.shape[0],moddt.shape[2],moddt.shape[3]))
    for e in range(moddt.shape[0]):
        for i in range(moddt.shape[2]):
            for j in range(moddt.shape[3]):
                xx = obsdt[:,i,j]
                yy = moddt[e,:,i,j]
                
                if np.isfinite(np.nanmax(xx)) == True and np.isfinite(np.nanmax(yy)) == True:
                    na = np.logical_or(np.isnan(xx),np.isnan(yy))
                    corr[e,i,j],pvals[e,i,j]=sts.pearsonr(xx[~na],yy[~na])
                else:
                    corr[e,i,j] = np.nan
                    pvals[e,i,j] = np.nan
                    
    ### Calculate ensemble mean
    meancorr = np.nanmean(corr,axis=0)
     
    #######################################################################
    #######################################################################
    #######################################################################
    ### Plot subplot of mean climate
    if variq == 'T2M':
        limit = np.arange(-1,1.01,0.1)
        barlim = np.round(np.arange(-1,1.1,0.5),2)
        cmap = cmocean.cm.balance
        label = r'\textbf{Anomaly Correlation Coefficient}'
        
    fig = plt.figure()
    var = meancorr
    
    ax1 = plt.subplot(111)
       
    m = Basemap(projection='ortho',lon_0=250,lat_0=40,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='dimgrey',linewidth=1)
    m.drawstates(color='dimgrey',linewidth=0.5)
    m.drawcountries(color='dimgrey',linewidth=0.5)
    
    circle = m.drawmapboundary(fill_color='white',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2,lat2,var,limit,extend='neither',latlon=True)
    
    cs1.set_cmap(cmocean.cm.balance)
    
    ## Box 1
    la1 = 43
    la2 = 60
    lo1 = 240
    lo2 = 295
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='blue', linewidth=2, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='blue', linewidth=2, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=2,color='blue',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=2,color='blue',zorder=4)
        
    fig.suptitle(r'\textbf{%s %s-%s -- SPEAR-AMIP}' % (yearmin,yearmax,slicemonth),fontsize=18,color='k')
    
    ###############################################################################
    cbar_ax1 = fig.add_axes([0.07,0.09,0.2,0.035])                
    cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                        extend='neither',extendfrac=0.07,drawedges=False)
    cbar1.set_label(label,fontsize=9,color='k',labelpad=1.4)  
    cbar1.set_ticks(barlim)
    cbar1.set_ticklabels(list(map(str,barlim)))
    cbar1.ax.tick_params(axis='x', size=.001,labelsize=7,pad=2)
    cbar1.outline.set_edgecolor('dimgrey')
    
    plt.tight_layout()
    if detrend == False:
        plt.savefig(directoryfigure + 'ACC_%s_AMIPold_%s_%s-%s.png'% (slicemonth,variq,yearmin,yearmax),dpi=300)
    else: 
        plt.savefig(directoryfigure + 'ACC_%s_AMIPold_%s_%s-%s_detrend.png'% (slicemonth,variq,yearmin,yearmax),dpi=300)
