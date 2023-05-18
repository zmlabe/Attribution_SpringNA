"""
Plot composites of AMIP data onto PNA/PCHZ1000 patterns

Author    : Zachary M. Labe
Date      : 12 July 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts
import read_FACTS_AMIPS as AM
import calc_DetrendData as DTT
import cmocean

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/validationSPEAR/RegressionPatternsAMIPS/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/' 
directorydataoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsAMIPS/'

mon_index = 'FMA'
mon_t2m = 'FMA'
detrend = True
variq = 'T2M'
yearsall = np.arange(1979,2021+1,1)
yearmin = 1979
yearmax = 2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]
yearAllData = years

indicesdir = ['PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
indices = ['PNAmodified','NAOmodified','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
indicesvari = ['Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500']
models = ['ECHAM5','ESRL-CAM5','spear']
modelname = ['echam','cam','spear']
scenario = ['amip_obs_rf']
level = 'surface'

### Sort functions
def calc_anomalies(years,data):
    """ 
    Calculate anomalies
    """
    
    ### Baseline - 1981-2010
    if data.ndim == 3:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[yearqold,:,:],axis=0)
        anoms = data - climold
    elif data.ndim == 4:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:]
    elif data.ndim == 5:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:,:]
    elif data.ndim == 6:
        yearqold = np.where((years >= 1981) & (years <= 2010))[0]
        climold = np.nanmean(data[:,yearqold,:,:],axis=1)
        anoms = data - climold[:,np.newaxis,:,:,:,:]
    
    return anoms

### Read in new indices
for mm in range(len(models)):
    ### Read in AMIP indices
    directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/Teleconnections/'
    amip = np.load(directoryoutput + 'ClimateIndices_ValidationSPEAR_%s_1979-2019.npz' % mon_index)
    amipi = amip['%s' % modelname[mm]][:,:]
    
    PNAq = amipi[0,:,:]
    AOq = amipi[1,:,:]

    ### Read in AMIP data
    if models[mm] == 'spear':
        lat1q,lon1q,lev1q,varq,yearsq = AM.read_FACTS_Experi('spear',models[mm],
                                                     variq,mon_t2m,5,
                                                     'nan',level)
    else:
        lat1q,lon1q,lev1q,varq,yearsq = AM.read_FACTS_Experi(scenario[0],models[mm],
                                                     variq,mon_t2m,5,
                                                     'nan',level)
    
    yearAllDataq = np.where((yearsq >= yearAllData.min()) & (yearsq <= yearAllData.max()))[0]
    varq = varq[:,yearAllDataq,:,:]
    print(yearAllDataq)
    
    ### Calculate regional mean
    anom = calc_anomalies(years,varq)
    
    ### Detrend t2m data for regressions
    if detrend == True:
        amip_t2m_dt = DTT.detrendData(anom,'surface','yearly')
    else:
        amip_t2m_dt = anom
        
    ensall = []
    for ee in range(len(amip_t2m_dt)):
        ensallq = []
        for yr in range(years.shape[0]):
            if ((PNAq[ee,yr] < 0.2) and (AOq[ee,yr] > 0.2)):
                ensallq.append(amip_t2m_dt[ee,yr,:,:])
        ensall.append(np.nanmean(np.asarray(ensallq),axis=0))
    
    ### Calculate ensemble mean
    meanAll = np.nanmean(np.asarray(ensall),axis=0)

    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Plot regression Pattern
    limit = np.arange(-1,1.01,0.05)
    barlim = np.round(np.arange(-1,2,1),2)
    
    fig = plt.figure(figsize=(6,6))
    label = r'\textbf{Composite of T2M [$^{\circ}$C]}'
    
    ax = plt.subplot(111)
    
    var = meanAll
    
    m = Basemap(projection='ortho',lon_0=265,lat_0=70,resolution='l',area_thresh=10000)
    m.drawcoastlines(color='dimgrey',linewidth=1)
    m.drawstates(color='dimgrey',linewidth=0.5)
    m.drawcountries(color='dimgrey',linewidth=0.5)
      
    if models[mm] != 'spear':
        var, lons_cyclic = addcyclic(var, lon1q)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat1q)
        x, y = m(lon2d, lat2d)
        
        circle = m.drawmapboundary(fill_color='white',color='dimgray',
                          linewidth=0.7)
        circle.set_clip_on(False)
        
        cs1 = m.contourf(x,y,var,limit,extend='both')
        
        cs1.set_cmap(cmocean.cm.balance)
    else:
        lon2,lat2 = np.meshgrid(lon1q,lat1q)
        
        circle = m.drawmapboundary(fill_color='white',color='dimgray',
                          linewidth=0.7)
        circle.set_clip_on(False)
        
        cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
        
        cs1.set_cmap(cmocean.cm.balance)
    
    ## Box 2
    la1 = 43
    la2 = 60
    lo1 = 240
    lo2 = 295
    lonsslice = np.linspace(lo1,lo2,lo2-lo1+1)
    latsslice = np.ones(len(lonsslice))*la2
    m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
    latsslice = np.ones(len(lonsslice))*la1
    m.plot(lonsslice, latsslice, color='b', linewidth=1.5, latlon=True,zorder=4)
    m.drawgreatcircle(lo1, la1, lo1, la2,linewidth=1.5,color='b',zorder=4)
    m.drawgreatcircle(lo2, la2, lo2, la1,linewidth=1.5,color='b',zorder=4)
    
    plt.title(r'\textbf{Composite Positive_negPNA-posNAO(%s) for %s(%s) using %s.png}' % (mon_index,variq,mon_t2m,models[mm]),fontsize=8,color='k')
        
    cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
    cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
    cbar1.set_ticks(barlim)
    cbar1.set_ticklabels(list(map(str,barlim)))
    cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
    cbar1.outline.set_edgecolor('dimgrey')
    plt.tight_layout()
        
    plt.savefig(directoryfigure + 'CompositeAMIP_pattern_negPNA-posNAO-%s_%s-%s_%s.png' % (mon_index,variq,mon_t2m,models[mm]),dpi=300)
