"""
Plot composites of different NAO/PNA combinations during period of interest

Author    : Zachary M. Labe
Date      : 30 June 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/work/Zachary.Labe/Data/ClimateIndices/'

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]

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
    
    return anoms

### Parameters
varplot = ['T2M','SST','SLP','Z500','Z100','Z30']
for v in range(len(varplot)):
    variqM = varplot[v]
    sliceperiod = 'April'
    monq = sliceperiod
    yearsall = np.arange(1979,2021+1,1)
    sliceshape = 3
    slicenan = 'nan'
    addclimo = True
    detrend = True
    yearmin = 1979
    yearmax = 2019
    
    ### Read data
    latobs,lonobs,varn = ERA.read_ERA5_monthly1x1(variqM,directorydata,sliceperiod,
                                                  yearsall,sliceshape,addclimo,slicenan)
    lon2,lat2 = np.meshgrid(lonobs,latobs)
    
    ### Read only 1979-2019
    yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
    years = yearsall[yearq]
    var = varn[yearq,:,:]
    
    ### Calculate anomalies
    anoms = calc_anomalies(years,var)
    
    ### Detrend data
    if detrend == True:
        vardt = DT.detrendDataR(np.asarray(anoms),'surface','monthly')
    else:
        vardt = anoms
    
    ### Read in teleconnections
    NAOd = np.genfromtxt(directorydata2 + 'NAO/NAO_CPC_1950-2021.txt',unpack=True)
    years_NAO = NAOd[0,:]
    NAO_index = NAOd[1:,yearq]
    
    PNAd = np.genfromtxt(directorydata2 + 'PNA/PNA_CPC_1950-2021.txt',unpack=True)
    years_PNA = PNAd[0,:]
    PNA_index = PNAd[1:,yearq]
    
    ### Slice months
    if monq == 'FMA':
        NAO = np.nanmean(NAO_index[1:4,:],axis=0)
        PNA = np.nanmean(PNA_index[1:4,:],axis=0)
    elif monq == 'FM':
        NAO = np.nanmean(NAO_index[1:3,:],axis=0)
        PNA = np.nanmean(PNA_index[1:3,:],axis=0)
    elif monq == 'March':
        NAO = np.nanmean(NAO_index[2:3,:],axis=0)
        PNA = np.nanmean(PNA_index[2:3,:],axis=0)
    elif monq == 'April':
        NAO = np.nanmean(NAO_index[3:4,:],axis=0)
        PNA = np.nanmean(PNA_index[3:4,:],axis=0)
    elif monq == 'none':
        NAO = NAO_index
        PNA = PNA_index
    else:
        print(ValueError('wrong months selected!'))
        sys.exit()
        
    posPNAposNAOq = np.where((PNA > 0) & (NAO > 0))[0]
    posPNAnegNAOq = np.where((PNA > 0) & (NAO < 0))[0]
    negPNAnegNAOq = np.where((PNA < 0) & (NAO < 0))[0]
    negPNAposNAOq = np.where((PNA < 0) & (NAO > 0))[0]
    # if len(posPNAposNAO) + len(posPNAnegNAO) + len(negPNAnegNAO) + len(negPNAposNAO) != years.shape[0]:
    #     print('SOMETHING IS WRONG WITH COMPOSITES!')
    #     sys.exit()
    
    ### Calculate mean composites
    posPNAposNAO = np.nanmean(vardt[posPNAposNAOq,:,:],axis=0)
    posPNAnegNAO = np.nanmean(vardt[posPNAnegNAOq,:,:],axis=0)
    negPNAnegNAO = np.nanmean(vardt[negPNAnegNAOq,:,:],axis=0)
    negPNAposNAO = np.nanmean(vardt[negPNAposNAOq,:,:],axis=0)
    
    ### Assemble for plotting
    nsample = [len(posPNAposNAOq),len(posPNAnegNAOq),len(negPNAposNAOq),len(negPNAnegNAOq)]
    allcomposites = [posPNAposNAO,posPNAnegNAO,negPNAposNAO,negPNAnegNAO]
    comptitles = ['+PNA,+NAO','+PNA,-NAO','-PNA,+NAO','-PNA,-NAO']
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Plot composite maps of different NAO/PNA combinations
    if variqM == 'T2M':
        limit = np.arange(-4,4.01,0.05)
        barlim = np.round(np.arange(-4,5,2),2)
        label = r'\textbf{%s Composites [$^{\circ}$C]}' % variqM
    elif variqM == 'SST':
        limit = np.arange(-1,1.01,0.01)
        barlim = np.round(np.arange(-1,2,1),2)
        label = r'\textbf{%s Composites [$^{\circ}$C]}' % variqM
    elif variqM == 'SLP':
        limit = np.arange(-4,4.01,0.1)
        barlim = np.round(np.arange(-4,5,2),2)
        label = r'\textbf{%s Composites [hPa]}' % variqM
    elif variqM == 'Z500':
        limit = np.arange(-40,40.01,5)
        barlim = np.round(np.arange(-40,50,20),2)
        label = r'\textbf{%s Composites [m]}' % variqM
    elif variqM == 'Z100':
        limit = np.arange(-200,200.01,20)
        barlim = np.round(np.arange(-200,300,100),2)
        label = r'\textbf{%s Composites [m]}' % variqM
    elif variqM == 'Z30':
        limit = np.arange(-200,200.01,20)
        barlim = np.round(np.arange(-200,300,100),2)
        label = r'\textbf{%s Composites [m]}' % variqM
    
    fig = plt.figure(figsize=(6,6))
    for i in range(len(allcomposites)):
        ax = plt.subplot(2,2,i+1)
        
        var = allcomposites[i]
        
        m = Basemap(projection='ortho',lon_0=265,lat_0=70,resolution='l',area_thresh=10000)
        m.drawcoastlines(color='dimgrey',linewidth=1)
        m.drawstates(color='dimgrey',linewidth=0.5)
        m.drawcountries(color='dimgrey',linewidth=0.5)
            
        var, lons_cyclic = addcyclic(var, lonobs)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, latobs)
        x, y = m(lon2d, lat2d)
        
        circle = m.drawmapboundary(fill_color='white',color='dimgray',
                          linewidth=0.7)
        circle.set_clip_on(False)
        
        cs1 = m.contourf(x,y,var,limit,extend='both')
        cs1.set_cmap(cmocean.cm.balance)
        
        if variqM == 'SST':
            m.fillcontinents(color='dimgrey')
            m.drawcoastlines(color='darkgrey',linewidth=0.5)
        
        plt.title(r'\textbf{%s [n=%s]}' % (comptitles[i],nsample[i]),fontsize=11,color='dimgrey')
        ax.annotate(r'\textbf{[%s]}' % (letters[i]),xy=(0,0),xytext=(0.07,0.90),
                  textcoords='axes fraction',color='k',fontsize=9,
                  rotation=0,ha='center',va='center')
        
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
        
    fig.suptitle(r'\textbf{%s-%s for %s}' % (yearmin,yearmax,sliceperiod),fontsize=15,color='k')
            
    cbar_ax1 = fig.add_axes([0.40,0.06,0.2,0.025])                
    cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
    cbar1.set_ticks(barlim)
    cbar1.set_ticklabels(list(map(str,barlim)))
    cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
    cbar1.outline.set_edgecolor('dimgrey')
    plt.tight_layout()
      
    if detrend == True:  
        plt.savefig(directoryfigure + 'CompositeMaps_%s-%s_NAO-PNA-%s_%s-%s_detrend.png' % (yearmin,yearmax,monq,variqM,sliceperiod),dpi=300)
    else:
        plt.savefig(directoryfigure + 'CompositeMaps_%s-%s_NAO-PNA-%s_%s-%s.png' % (yearmin,yearmax,monq,variqM,sliceperiod),dpi=300)

