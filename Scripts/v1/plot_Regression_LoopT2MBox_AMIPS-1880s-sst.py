"""
Plot regression patterns for amips using the T2M_BoxNA (spear - global sst maps)

Author    : Zachary M. Labe
Date      : 22 July 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/AMIPS/RegressionPatternsAMIPS/T2M_BoxNA/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 
directorydataoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsAMIPS/T2M_BoxNA/'

variregress = ['SST']
monthsall = ['January','February','March','April','JFM','FM','FMA']
models = ['SPEAR']
modelname = ['echam','cam','spear']
scenario = ['amip_obs_rf']
level = 'surface'

for mo in range(len(monthsall)):
    for moo in range(len(monthsall)):
        mon_index = monthsall[mo]
        mon_t2m = monthsall[moo]
        detrend = True
        
        for iii in range(len(variregress)):
            variq = variregress[iii]
            yearsall = np.arange(1979,2021+1,1)
            yearmin = 1979
            yearmax = 2019
            yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
            years = yearsall[yearq]
            yearAllData = years

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
            
            def regressData(x,y,runnamesm):    
                print('\n>>> Using regressData function!')
                
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
            
            ### Read in new indices
            for mm in range(len(models)):
                ### Read in AMIP indices
                directoryoutput = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/Teleconnections/'
                amip = np.load(directoryoutput + 'ClimateIndices-T2M_BoxNA_AMIPSPEAR_%s_1979-2019.npz' % mon_index)
                amipi = amip['%s' % modelname[mm]][0,:] # only one index in the file
            
                ### Read in AMIP data
                lat1q,lon1q,lev1q,varq,yearsq = AM.read_FACTS_Experi('spear_obs_rf',models[mm],
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
                    
                ### Calculate regression
                ensRegr = np.empty((len(anom),lat1q.shape[0],lon1q.shape[0]))
                ensPVal = np.empty((len(anom),lat1q.shape[0],lon1q.shape[0]))
                ensRVal = np.empty((len(anom),lat1q.shape[0],lon1q.shape[0]))
                for ee in range(len(anom)):
                    slope,intercept,rvalue,pvalue,stderr = regressData(amipi[ee],np.asarray(amip_t2m_dt[ee]),'obs')
                    ensRegr[ee,:,:] = slope
                    ensPVal[ee,:,:] = pvalue
                    ensRVal[ee,:,:] = rvalue
                    print('>>> Finished calculating ensemble-%s for regression of %s!' % (ee+1,models[mm]))
                
                ### Calculate ensemble means
                slopeM = np.nanmean(ensRegr,axis=0)
                pvalueM = np.nanmean(ensPVal,axis=0)
                
                ### Significant at 95% confidence level
                pvalueM[np.where(pvalueM >= 0.05)] = np.nan
                pvalueM[np.where(pvalueM < 0.05)] = 1.
                ###############################################################################
                ###############################################################################
                ###############################################################################
                ### Plot regression Pattern
                if variq == 'SST':
                    limit = np.arange(-0.50,0.501,0.01)
                    barlim = np.round(np.arange(-0.5,0.51,0.25),2)
                
                fig = plt.figure(figsize=(6,6))
                label = r'\textbf{Regression Coefficients}'
                
                ax = plt.subplot(111)
                
                var = slopeM
                pvar = pvalueM
                
                m = Basemap(projection='robin',lon_0=-180,resolution='l',area_thresh=10000)
                m.drawcoastlines(color='dimgrey',linewidth=1)
                m.drawstates(color='dimgrey',linewidth=0.5)
                m.drawcountries(color='dimgrey',linewidth=0.5)

                var,lon1qs = shiftgrid(180,var,lon1q,start=False)
                  

                lon2,lat2 = np.meshgrid(lon1qs,lat1q)
                
                circle = m.drawmapboundary(fill_color='white',color='dimgray',
                                  linewidth=0.7)
                circle.set_clip_on(False)
                
                cs1 = m.contourf(lon2,lat2,var,limit,extend='both',latlon=True)
                cs2 = m.contourf(lon2,lat2,pvar,colors='None',hatches=['....'],
                             linewidths=0.4,latlon=True)
                
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
                
                plt.title(r'\textbf{Regression Pattern for %s(%s) on %s(%s) using %s}' % (variq,mon_t2m,'T2M_BoxNA',mon_index,models[mm]),fontsize=8,color='k')
                    
                cbar_ax1 = fig.add_axes([0.05,0.06,0.2,0.025])                
                cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                                    extend='max',extendfrac=0.07,drawedges=False)
                cbar1.set_label(label,fontsize=7,color='dimgrey',labelpad=1.4)  
                cbar1.set_ticks(barlim)
                cbar1.set_ticklabels(list(map(str,barlim)))
                cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
                cbar1.outline.set_edgecolor('dimgrey')
                plt.tight_layout()
                    
                plt.savefig(directoryfigure + 'RegressionAMIP_T2M_BoxNApattern_%s-%s_%s-%s_%s_globalsMaps.png' % (variq,mon_t2m,'T2M_BoxNA',mon_index,models[mm]),dpi=300)
