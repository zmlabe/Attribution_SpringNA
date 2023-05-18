"""
Calculate regression maps for loops of indices onto T2MBox_NA

Author    : Zachary M. Labe
Date      : 17 November 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/ObsRegressions/T2M_BoxNA/' 
directorydata = '/work/Zachary.Labe/Data/' 
directorydata2 = '/work/Zachary.Labe/Data/ClimateIndices/'
directorydataoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/Obs/'

### Loop through indices
variregress = ['SND','Z500','SLP','SST','T2M','Z30','Z100','Z1000','U200','U700']

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
monthsall = ['January','February','March','April','May','JFM','FM','FMA']
for mo in range(len(monthsall)):
    for moo in range(len(monthsall)):
        monq = '%s' % monthsall[mo] # earlier
        
        sliceperiod = monthsall[moo] # later
        yearsall = np.arange(1979,2021+1,1)
        sliceshape = 3
        slicenan = 'nan'
        addclimo = True
        yearmin = 1979
        yearmax = 2020
    
        for iii in range(len(variregress)):
            variqM = variregress[iii]
        
            ### Read data
            latobs,lonobs,levobs,varn = ERA.read_ERA5_monthly1x1(variqM,directorydata,sliceperiod,
                                                          yearsall,sliceshape,addclimo,slicenan,'surface')
            lon2,lat2 = np.meshgrid(lonobs,latobs)
            
            ### Read only 1979-2020
            yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
            years = yearsall[yearq]
            var = np.asarray(varn[yearq,:,:])
            
            ### Calculate anomalies
            anoms = calc_anomalies(years,var)
            
            ### Detrend data
            vardt = DT.detrendDataR(anoms,'surface','monthly')
            
            ### Read in teleconnections
            T2Md = np.genfromtxt(directorydata2 + '%s/%s_%s_%s_%s-%s_detrended.txt' % ('T2M_BoxNA','T2M_BoxNA','T2M',monq,yearmin,yearmax),unpack=True)
            years_T2M = T2Md[1,:]
            yearclimoq = np.where((years>=1981) & (years<=2010))[0]
            years_T2Mq = years_T2M[yearclimoq]
            meanclimo = np.nanmean(years_T2Mq)
            stdclimo = np.nanstd(years_T2Mq)
            T2M = (years_T2M - meanclimo)/stdclimo
            
            ### Calculate regression
            slope,intercept,rvalue,pvalue,stderr = regressData(T2M,vardt,'obs')
            
            ### Significant at 95% confidence level
            pvalue[np.where(pvalue >= 0.05)] = np.nan
            pvalue[np.where(pvalue < 0.05)] = 1.
            
            ###############################################################################
            ###############################################################################
            ###############################################################################
            ### Plot subplot of different SAI analysis
            if variqM == 'T2M':
                limit = np.arange(-2,2.01,0.05)
                barlim = np.round(np.arange(-2,3,1),2)
            elif variqM == 'SST':
                limit = np.arange(-0.50,0.501,0.01)
                barlim = np.round(np.arange(-0.5,0.51,0.25),2)
            elif variqM == 'SLP':
                limit = np.arange(-5,5.01,0.05)
                barlim = np.round(np.arange(-5,6,1),2)
            elif variqM == 'SND':
                limit = np.arange(-0.04,0.041,0.001)
                barlim = np.round(np.arange(-0.04,0.05,0.04),2)
            elif any([variqM == 'Z30',variqM == 'Z100',variqM == 'Z1000',variqM == 'Z500']):
                limit = np.arange(-50,50.01,0.5)
                barlim = np.round(np.arange(-50,51,50),2)
            elif any([variqM == 'U200']):
                limit = np.arange(-5,5.01,0.05)
                barlim = np.round(np.arange(-5,6,1),2)
            elif any([variqM == 'U200']):
                limit = np.arange(-2,2.01,0.05)
                barlim = np.round(np.arange(-2,3,1),2)
            
            fig = plt.figure(figsize=(6,6))
            label = r'\textbf{Regression Coefficients}'
            
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
            
            plt.title(r'\textbf{%s-%s of %s(%s)-%s(%s)}' % (yearmin,yearmax,'T2M_BoxNA',monq,variqM,sliceperiod),fontsize=15,color='k')
                
            cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
            cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                                extend='max',extendfrac=0.07,drawedges=False)
            cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
            cbar1.set_ticks(barlim)
            cbar1.set_ticklabels(list(map(str,barlim)))
            cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
            cbar1.outline.set_edgecolor('dimgrey')
            plt.tight_layout()
        
            plt.savefig(directoryfigure + 'RegressionObs_T2M_BoxNApattern_%s-%s_%s-%s.png' % (variqM,sliceperiod,'T2M_BoxNA',monq),dpi=300)
            
            ### Save file
            np.savez(directorydataoutput + 'RegressionObs_T2M_BoxNApattern_%s-%s_%s-%s.npz' % (variqM,sliceperiod,'T2M_BoxNA',monq),
                     regress=slope,pval=pvalue,r=rvalue,lat=latobs,lon=lonobs)
