"""
Evaluate relationships in stratosphere for CMIP6 models

Author    : Zachary M. Labe
Date      : 12 August 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import cmocean
import numpy as np
import calc_Utilities as UT
import sys
import read_CMIP6 as CM
import calc_DetrendData as DT
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/CMIP6/' 
directorydata = '/work/Zachary.Labe/Data/'

### Functions
def readInVar(variq,GCM,months,sliceshape,periodyrs,slinenan,level):
    cmipmodels = []
    for i in range(len(GCM)):
        lat1,lon1,lev,var = CM.read_CMIP6_monthly(variq,GCM[i],months,sliceshape,
                                                  periodyrs,slicenan,level)
        cmipmodels.append(var)
    cmipmodels = np.asarray(cmipmodels).squeeze()
    lon2,lat2 = np.meshgrid(lon1,lat1)
    return lat1,lon1,lat2,lon2,lev,cmipmodels

def calc_regionalAve(regionn,lat1,lon1,data):
    print('\n>>> CALCULATING REGIONAL AVERAGE!')
    
    ### Pick region
    if regionn == 'CANMID':
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]      
    elif regionn == 'PCH':
        la1 = 65
        la2 = 90
        lo1 = 0
        lo2 = 360
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]           
    else:
        print(ValueError('WRONG REGION SELECTED!'))
        sys.exit()
        
    ### Slice region
    if data.ndim == 3:
        meanlat = data[:,lat1q,:]
        meanbox = meanlat[:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    elif data.ndim == 4:
        meanlat = data[:,:,lat1q,:]
        meanbox = meanlat[:,:,:,lon1q]
        lon1a = lon1[lon1q]
        lat1a = lat1[lat1q]
        lon2q,lat2q = np.meshgrid(lon1a,lat1a)
        
    ### Calculate timeseries
    print('Shape is --->',(meanbox.shape))
    mean = UT.calc_weightedAve(meanbox,lat2q)

    print('>>> DONE REGIONAL AVERAGE!')
    return mean,lat1a,lon1a,lat2q,lon2q

def regressData(x,y,runnamesm):    
    print('\n>>> Using regressData function! \n')
    
    if y.ndim == 5: # 5D array
        slope = np.empty((y.shape[0],y.shape[1],y.shape[3],y.shape[4]))
        intercept = np.empty((y.shape[0],y.shape[1],y.shape[3],y.shape[4]))
        rvalue = np.empty((y.shape[0],y.shape[1],y.shape[3],y.shape[4]))
        pvalue = np.empty((y.shape[0],y.shape[1],y.shape[3],y.shape[4]))
        stderr = np.empty((y.shape[0],y.shape[1],y.shape[3],y.shape[4]))
        for model in range(y.shape[0]):
            print('Completed: Regression for %s!' % runnamesm[model])
            for ens in range(y.shape[1]):
                for i in range(y.shape[3]):
                    for j in range(y.shape[4]):
                        ### 1D time series for regression
                        xx = x
                        yy = y[model,ens,:,i,j]
                        
                        ### Mask data for nans
                        mask = ~np.isnan(xx) & ~np.isnan(yy)
                        varx = xx[mask]
                        vary = yy[mask]
                        
                        ### Calculate regressions
                        slope[model,ens,i,j],intercept[model,ens,i,j], \
                        rvalue[model,ens,i,j],pvalue[model,ens,i,j], \
                        stderr[model,ens,i,j] = sts.linregress(varx,vary)
                        
    if y.ndim == 4: # 4D array
        slope = np.empty((y.shape[0],y.shape[2],y.shape[3]))
        intercept = np.empty((y.shape[0],y.shape[2],y.shape[3],))
        rvalue = np.empty((y.shape[0],y.shape[2],y.shape[3]))
        pvalue = np.empty((y.shape[0],y.shape[2],y.shape[3],))
        stderr = np.empty((y.shape[0],y.shape[2],y.shape[3]))
        for model in range(y.shape[0]):
            print('Completed: Regression for %s!' % runnamesm[model])
            for i in range(y.shape[2]):
                for j in range(y.shape[3]):
                    ### 1D time series for regression
                    xx = x
                    yy = y[model,:,i,j]
                    
                    ### Mask data for nans
                    mask = ~np.isnan(xx) & ~np.isnan(yy)
                    varx = xx[mask]
                    vary = yy[mask]
                        
                    ### Calculate regressions
                    slope[model,i,j],intercept[model,i,j], \
                    rvalue[model,i,j],pvalue[model,i,j], \
                    stderr[model,i,j] = sts.linregress(varx,vary)
                    
    elif y.ndim == 3: #3D array
        slope = np.empty((y.shape[1],y.shape[2]))
        intercept = np.empty((y.shape[1],y.shape[2]))
        rvalue = np.empty((y.shape[1],y.shape[2]))
        pvalue = np.empty((y.shape[1],y.shape[2]))
        stderr = np.empty((y.shape[1],y.shape[2]))
        for i in range(y.shape[1]):
            for j in range(y.shape[2]):
                ### 1D time series for regression
                xx = x
                yy = y[:,i,j]
                
                ### Mask data for nans
                mask = ~np.isnan(xx) & ~np.isnan(yy)
                varx = xx[mask]
                vary = yy[mask]
                        
                ### Calculate regressions
                if np.isfinite(np.nanmean(yy)):
                    slope[i,j],intercept[i,j],rvalue[i,j], \
                        pvalue[i,j],stderr[i,j] = sts.linregress(varx,vary)
                else:
                    slope[i,j] = np.nan
                    intercept[i,j] = np.nan
                    rvalue[i,j] = np.nan
                    pvalue[i,j] = np.nan
                    stderr[i,j] = np.nan
                        
    print('>>> Completed: Finished regressData function!')
    return slope,intercept,rvalue,pvalue,stderr

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
modelsall = ['CESM2_WACCM6','MIROC6']
mont_t = 'FMA'
mont_ind = 'FMA'
yearsall = np.arange(1979,2022+1,1)
periodyrs = '1979-2022'
sliceshape = 4
slicenan = 'nan'
level = 'surface'
variq_t = 'T2M'
variq_ind = 'Z30'

### Read in models
lat1,lon1,lat2,lon2,lev,tas = readInVar(variq_t,modelsall,mont_t,sliceshape,
                                        periodyrs,slicenan,level)
lat1,lon1,lat2,lon2,lev,z100 = readInVar(variq_ind,modelsall,mont_ind,sliceshape,
                                        periodyrs,slicenan,level)

### Calculate anomalies
tas_anom = calc_anomalies(yearsall,tas)
z100_anom = calc_anomalies(yearsall,z100)

### Detrend data at each grid point
tas_dt = DT.detrendData(tas_anom,'surface',yearsall)
z100_dt = DT.detrendData(z100_anom,'surface',yearsall)

### Calculate regional averages
tas_can,lat1a,lon1a,lat2a,lon2a = calc_regionalAve('CANMID',lat1,lon1,tas_dt)
z100_pch,lat1a,lon1a,lat2a,lon2a = calc_regionalAve('PCH',lat1,lon1,z100_dt)

### Standardize index
tas_canz = sts.zscore(tas_can,axis=1)
z100_pchz = sts.zscore(z100_pch,axis=1)

### Calculate correlations
corr_waccm = sts.pearsonr(tas_canz[0],z100_pchz[0])
corr_miroc = sts.pearsonr(tas_canz[1],z100_pchz[1])

### Calculate regression
slope = np.empty((len(modelsall),lat1.shape[0],lon1.shape[0]))
intercept = np.empty((len(modelsall),lat1.shape[0],lon1.shape[0]))
rvalue = np.empty((len(modelsall),lat1.shape[0],lon1.shape[0]))
pvalue = np.empty((len(modelsall),lat1.shape[0],lon1.shape[0]))
stderr = np.empty((len(modelsall),lat1.shape[0],lon1.shape[0]))
for m in range(len(modelsall)):
    slopeq,interceptq,rvalueq,pvalueq,stderrq = regressData(z100_pchz[m],tas_dt[m],modelsall[m])

    ### Significant at 95% confidence level
    pvalueq[np.where(pvalueq >= 0.05)] = np.nan
    pvalueq[np.where(pvalueq < 0.05)] = 1.
    
    slope[m,:,:] = slopeq
    intercept [m,:,:] = interceptq
    rvalue[m,:,:] = rvalueq
    pvalue[m,:,:] = pvalueq
    stderr[m,:,:] = stderrq
    
###############################################################################
###############################################################################
###############################################################################
### Plot subplots of regressions
limit = np.arange(-1,1.01,0.05)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(8,4))
label = r'\textbf{Regression Coefficients [$^{\circ}$C/m]}'

ax = plt.subplot(121)

var = slope[0]
pvar = pvalue[0]

m = Basemap(projection='ortho',lon_0=265,lat_0=70,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)
    
pvar,lons_cyclic = addcyclic(pvar, lon1)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,var,limit,extend='both')
cs2 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
             linewidths=0.4)

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

plt.title(r'\textbf{%s for %s-%s -- %s-%s}' % (modelsall[0],variq_ind,mont_ind,variq_t,mont_t),
          fontsize=11,color='dimgrey')

###############################################################################
###############################################################################
###############################################################################
ax = plt.subplot(122)

var = slope[1]
pvar = pvalue[1]

m = Basemap(projection='ortho',lon_0=265,lat_0=70,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)
    
pvar,lons_cyclic = addcyclic(pvar, lon1)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
var, lons_cyclic = addcyclic(var, lon1)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

circle = m.drawmapboundary(fill_color='white',color='dimgray',
                  linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,var,limit,extend='both')
cs2 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
             linewidths=0.4)

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

plt.title(r'\textbf{%s for %s-%s -- %s-%s}' % (modelsall[1],variq_ind,mont_ind,variq_t,mont_t),
          fontsize=11,color='dimgrey')
    
cbar_ax1 = fig.add_axes([0.401,0.08,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()

plt.savefig(directoryfigure + 'CMIP6_%s-%s_%s-%s.png' % (variq_ind,mont_ind,variq_t,mont_t),dpi=300)
