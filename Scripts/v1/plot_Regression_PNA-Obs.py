"""
Calculate regression maps of PNA onto T2M maps

Author    : Zachary M. Labe
Date      : 29 June 2022
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
                slope[i,j],intercept[i,j],rvalue[i,j], \
                pvalue[i,j],stderr[i,j] = sts.linregress(varx,vary)
                        
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
variqM = 'T2M'
sliceperiod = 'FMA' # later
monq = 'FMA' # earlier
yearsall = np.arange(1979,2021+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
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
vardt = DT.detrendDataR(anoms,'surface','monthly')

### Read in teleconnections
aod = np.genfromtxt(directorydata2 + 'AO/AO_CPC_1950-2021.txt',unpack=True)
years_ao = aod[0,:]
ao_index = aod[1:,yearq]

naod = np.genfromtxt(directorydata2 + 'NAO/NAO_CPC_1950-2021.txt',unpack=True)
years_nao = naod[0,:]
nao_index = naod[1:,yearq]

pnad = np.genfromtxt(directorydata2 + 'PNA/PNA_CPC_1950-2021.txt',unpack=True)
years_pna = pnad[0,:]
pna_index = pnad[1:,yearq]

### Slice months
if monq == 'FMA':
    AO = np.nanmean(ao_index[1:4,:],axis=0)
    NAO = np.nanmean(nao_index[1:4,:],axis=0)
    PNA = np.nanmean(pna_index[1:4,:],axis=0)
elif monq == 'FM':
    AO = np.nanmean(ao_index[1:3,:],axis=0)
    NAO = np.nanmean(nao_index[1:3,:],axis=0)
    PNA = np.nanmean(pna_index[1:3,:],axis=0)
elif monq == 'March':
    AO = np.nanmean(ao_index[2:3,:],axis=0)
    PNA = np.nanmean(pna_index[2:3,:],axis=0)
    NAO = np.nanmean(nao_index[2:3,:],axis=0)
elif monq == 'April':
    AO = np.nanmean(ao_index[3:4,:],axis=0)
    PNA = np.nanmean(pna_index[3:4,:],axis=0)
    NAO = np.nanmean(nao_index[3:4,:],axis=0)
elif monq == 'none':
    AO = ao_index
    PNA = pna_index
    NAO = nao_index
else:
    print(ValueError('wrong months selected!'))
    sys.exit()
    
### Calculate regression
slope,intercept,rvalue,pvalue,stderr = regressData(PNA,vardt,'obs')

### Significant at 95% confidence level
pvalue[np.where(pvalue >= 0.05)] = np.nan
pvalue[np.where(pvalue < 0.05)] = 1.

###############################################################################
###############################################################################
###############################################################################
### Plot subplot of different SAI analysis
limit = np.arange(-1,1.01,0.05)
barlim = np.round(np.arange(-1,2,1),2)

fig = plt.figure(figsize=(6,6))
label = r'\textbf{Regression Coefficients [$^{\circ}$C/m]}'

ax = plt.subplot(111)

var = slope
pvar = pvalue

m = Basemap(projection='ortho',lon_0=265,lat_0=70,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)
    
pvar,lons_cyclic = addcyclic(pvar, lonobs)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
var, lons_cyclic = addcyclic(var, lonobs)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, latobs)
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

plt.title(r'\textbf{%s-%s of PNAi(%s)-T2M(%s)}' % (yearmin,yearmax,monq,sliceperiod),fontsize=15,color='k')
    
cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')
plt.tight_layout()
    
plt.savefig(directoryfigure + 'RegressionMap_%s-%s_PNAi-%s_T2M-%s.png' % (yearmin,yearmax,monq,sliceperiod),dpi=300)

