"""
Calculate regression maps for loops of indices

Author    : Zachary M. Labe
Date      : 17 July 2022
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
from scipy.interpolate import griddata as g

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/ObsRegressions/ModelComparison/' 
directorydataObs = '/work/Zachary.Labe/Data/ClimateIndices/'
directorydataObsT2Mbox = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/'

### Parameters
yearsall = np.arange(1979,2021+1,1)
yearmin = 1979
yearmax = 2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

### Loop through indices
# indicesdir = ['PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI','NINO34','WPPRECT','NECPSST']
# indices = ['PNAmodified','NAOmodified','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI','NINO34','WPPRECT','NECPSST']
# indicesvari = ['Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500','SST','P','SST']
indicesdir = ['PNA','NAO','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
indices = ['PNAmodified','NAOmodified','NPI','PCHT2M','PCHZ30','PCHZ100','PCHZ1000','SHI','UBI']
indicesvari = ['Z500','SLP','SLP','T2M','Z30','Z100','Z1000','SLP','Z500']
monthsloop = ['FMA']
t2mmonth = 'FMA'
model = 'ESRL-CAM5'
region = 'NA'

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

def regrid(lat11,lon11,lat21,lon21,var):
    """
    Interpolated on selected grid. Reads ERA5 in as 4d with 
    [year,month,lat,lon]
    """
    
    lon1,lat1 = np.meshgrid(lon11,lat11)
    lon2,lat2 = np.meshgrid(lon21,lat21)
    
    varn_re = np.reshape(var,((lat1.shape[0]*lon1.shape[1])))   
    
    print('Completed: Start regridding process:')
    z = g((np.ravel(lat1),np.ravel(lon1)),varn_re,(lat2,lon2),method='linear')
    print('Completed: Regridding---')
    return z

### Function to calculate pattern correlations over different regions
def patternCorrelationRegion(dataobs,datamodel,lat1,lon1,region):
    if region != 'global':
        if region == 'focusarea':
            latmin = 43
            latmax = 60
            lonmin = 240
            lonmax = 295
        elif region == 'arctic':
            latmin = 65
            latmax = 90
            lonmin = 0
            lonmax = 360
        elif region == 'NH':
            latmin = 0
            latmax = 90
            lonmin = 0
            lonmax = 360
        elif region == 'NHextra':
            latmin = 30
            latmax = 90
            lonmin = 0
            lonmax = 360
        elif region == 'NA':
            latmin = 10
            latmax = 80
            lonmin = 180
            lonmax = 300
        elif region == 'US':
            latmin = 24
            latmax = 60
            lonmin = 235
            lonmax = 290
        
        latq = np.where((lat1 > latmin) & (lat1 < latmax))[0]
        latnew = lat1[latq]
        dataobs_new1 = dataobs[latq,:]
        datamodels_new1 = datamodel[latq,:]
    
        if np.max(lon1) < 200:
            print(ValueError('Longitude grid is not correct!'))
            sys.exit()
        lonq = np.where((lon1 > lonmin) & (lon1 < lonmax))[0]
        lonnew = lon1[lonq]
        dataobs_new2 = dataobs_new1[:,lonq]
        datamodels_new2 = datamodels_new1[:,lonq]

        ### Prepare for pattern correlations
        dataobs_new = dataobs_new2
        datamodels_new = datamodels_new2
    
    else:
        latnew = lat1
        lonnew = lon1
        dataobs_new = dataobs
        datamodels_new = datamodel
    
    ### Calculate pattern correlation
    corr = UT.calc_spatialCorr(dataobs_new,datamodels_new,latnew,lonnew,'yesnan')
    
    return corr,latnew,lonnew

for mmm in range(len(monthsloop)):
    corriii = []
    for iii in range(len(indices)):
        regressobsq = np.load('/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsObs/RegressionObs_pattern_%s-%s_%s-%s.png.npz' % ('T2M',t2mmonth,indices[iii],monthsloop[mmm]))
        regressobs = regressobsq['regress']
        latobs = regressobsq['lat']
        lonobs = regressobsq['lon']
       
        if model == 'coupledSPEAR':
            if iii < 9:
                regressmodelq = np.load('/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsCoupledSPEAR/RegressionCoupledSPEAR_pattern_%s-%s_%s-%s_SPEAR-MED.npz' % ('T2M',t2mmonth,indices[iii],monthsloop[mmm]))
                regressmodel = regressmodelq['regress']
                latmodel = regressmodelq['lat']
                lonmodel = regressmodelq['lon']   
            else:
                regressmodelq = np.load('/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsCoupledSPEAR/RegressionCoupledSPEAR_Tropics_pattern_%s-%s_%s-%s_SPEAR-MED.npz' % ('T2M',t2mmonth,indices[iii],monthsloop[mmm]))
                regressmodel = regressmodelq['regress']
                latmodel = regressmodelq['lat']
                lonmodel = regressmodelq['lon']    
        elif model == 'coupledSPEAR-NOAER':
            if iii < 9:
                regressmodelq = np.load('/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsCoupledSPEAR-NOAER/RegressionCoupledSPEAR-NOAER_pattern_%s-%s_%s-%s_SPEAR-MED-NOAER.npz' % ('T2M',t2mmonth,indices[iii],monthsloop[mmm]))
                regressmodel = regressmodelq['regress']
                latmodel = regressmodelq['lat']
                lonmodel = regressmodelq['lon']   
            else:
                regressmodelq = np.load('/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsCoupledSPEAR-NOAER/RegressionCoupledSPEAR-NOAER_Tropics_pattern_%s-%s_%s-%s_SPEAR-MED-NOAER.npz' % ('T2M',t2mmonth,indices[iii],monthsloop[mmm]))
                regressmodel = regressmodelq['regress']
                latmodel = regressmodelq['lat']
                lonmodel = regressmodelq['lon']   
        elif any([model == 'ECHAM5',model=='spear',model=='ESRL-CAM5']):
            regressmodelq = np.load('/home/Zachary.Labe/Research/Attribution_SpringNA/Data/RegressionPatternsAMIPS/RegressionAMIP_pattern_%s-%s_%s-%s_%s.png.npz' % ('T2M',t2mmonth,indices[iii],monthsloop[mmm],model))
            regressmodel = regressmodelq['regress']
            latmodel = regressmodelq['lat']
            lonmodel = regressmodelq['lon']   
                
        ### Regrid to match models
        newobs = regrid(latobs,lonobs,latmodel,lonmodel,regressobs)                  

        corr = np.empty((len(regressmodel)))
        for e in range(len(regressmodel)):
            corr[e],latnew,lonnew = patternCorrelationRegion(newobs,regressmodel[e],latmodel,lonmodel,region)
            
        ### Save for plotting
        corriii.append(corr)
    
    ### Make box plots
    corall = np.asarray(corriii)
    
    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(111)
    
    plotdata = corall
    
    adjust_spines(ax, ['left', 'bottom'])
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_linewidth(2)
    ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
    ax.tick_params(axis="x",which="both",bottom = False,top=False,
                    labelbottom=True)
    
    ax.yaxis.grid(zorder=1,color='darkgrey',alpha=0.7,clip_on=False,linewidth=0.5)
    
    def set_box_color(bp, color):
        plt.setp(bp['boxes'],color=color)
        plt.setp(bp['whiskers'], color=color,linewidth=1.5)
        plt.setp(bp['caps'], color='w',alpha=0)
        plt.setp(bp['medians'], color='w',linewidth=1)
    
    positionsq = np.arange(len(corall))
    bpl = plt.boxplot(plotdata.transpose(),positions=positionsq,widths=0.6,
                      patch_artist=True,sym='')
    
    # Modify boxes
    cp= 'maroon'
    set_box_color(bpl,cp)
    plt.plot([], c=cp, label=r'\textbf{PATTERN CORRELATION}',clip_on=False)
        
    for i in range(len(plotdata)):
        y = plotdata[i]
        x = np.random.normal(positionsq[i], 0.04, size=len(y))
        plt.plot(x, y,color='teal', alpha=0.5,zorder=10,marker='.',linewidth=0,markersize=5,markeredgewidth=0,clip_on=False)
     
    plt.yticks(np.arange(-1,1.1,0.2),list(map(str,np.round(np.arange(-1,1.1,0.2),2))),
                fontsize=6) 
    plt.ylim([-1,1])
    
    plt.xticks(np.arange(len(corall)),indicesdir,fontsize=6) 
    plt.xlim([-1,len(corall)+1])
    plt.title(r'\textbf{%s-%s on INDEX-%s using %s}' % ('T2M',t2mmonth,monthsloop[mmm],model),color='k',fontsize=10)
    
    plt.ylabel(r'\textbf{Pattern Correlation -- %s}' % monthsloop[mmm],color='k',fontsize=7)

    plt.savefig(directoryfigure + 'patternCorrelations_RegressionPatterns_%s-%s_INDEX-%s_%s_%s.png' % ('T2M',t2mmonth,monthsloop[mmm],model,region),dpi=300)

