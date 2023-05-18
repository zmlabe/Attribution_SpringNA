"""
Plot differences of climate baselines 1981-2010 and 1991-2020 by month
Author    : Zachary M. Labe
Date      : 16 May 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import read_ERA5_monthlyBE as ERA
import scipy.stats as sts
import read_GISTEMP as GG

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/Users/zlabe/Documents/Research/Attribution_SpringNA/Figures/' 
directorydata = '/Users/zlabe//Data/ERA5/' 
directorydata2 = '/Users/zlabe//Data/GISTEMP/'

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'April'
years = np.arange(1950,2019+1,1)
years2 = np.arange(1900,2019+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True
length = [10,30]

fig = plt.figure(figsize=(8,4))
for l in range(len(length)):
    ### Read data
    lat1,lon1,data = ERA.read_ERA5_monthlyBE(variq,directorydata,sliceperiod,
                                    years,sliceshape,addclimo,
                                    slicenan)
    
    data = data[:-2]
    
    ### Baseline - 1981-2010
    yearqold = np.where((years >= 1951) & (years <= 1980))[0]
    climold = np.nanmean(data[yearqold,:,:],axis=0)
    anommean = data - climold
    
    ### Calculate region
    la1 = 43
    la2 = 60
    lo1 = 240
    lo2 = 295
    lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
    lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
    
    meanlat = anommean[:,lat1q,:]
    meanbox = meanlat[:,:,lon1q]
    lon1a = lon1[lon1q]
    lat1a = lat1[lat1q]
    lon2q,lat2q = np.meshgrid(lon1a,lat1a)
    
    ### Calculate timeseries
    mean = UT.calc_weightedAve(meanbox,lat2q)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Read data
    lat1c,lon1c,varc = GG.read_GISTEMP(directorydata2,sliceperiod,years2,sliceshape,addclimo,slicenan,False)
    
    varc = varc[:-2]
    
    ### Calculate region
    lat1qc = np.where((lat1c >= la1) & (lat1c <= la2))[0]
    lon1qc = np.where((lon1c >= lo1) & (lon1c <= lo2))[0]
    
    meanlatc = varc[:,lat1qc,:]
    meanboxc = meanlatc[:,:,lon1q]
    lon1ac = lon1c[lon1qc]
    lat1ac = lat1c[lat1qc]
    lon2qc,lat2qc = np.meshgrid(lon1ac,lat1ac)
    
    ### Calculate timeseries
    meanc = UT.calc_weightedAve(meanboxc,lat2qc)
    
    def calcTrends(yearsnew,trendlength,data):
        ### Calculate trend periods
        yearstrend = np.empty((len(yearsnew)-trendlength+1,trendlength))
        datatrend = np.empty((len(yearsnew)-trendlength+1,trendlength))
        for hi in range(len(yearsnew)-(trendlength-1)):
            yearstrend[hi,:] = np.arange(yearsnew[hi],yearsnew[hi]+trendlength,1)
            datatrend[hi,:] = data[hi:hi+trendlength]
            
        print(yearstrend)
        
        ### Calculate trend lines    
        linetrend = np.empty((len(yearsnew)-trendlength+1,2))
        for hi in range(len(yearsnew)-trendlength+1):         
            linetrend[hi,:] = np.polyfit(yearstrend[hi],datatrend[hi],1)
            
        slope = linetrend[:,0]
        return slope
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Calculate slopes
    linetrend = calcTrends(years,length[l],mean)*10
    linetrendc = calcTrends(years2,length[l],meanc)*10
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
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
            
    ##############################################################################
    ##############################################################################
    ##############################################################################
    ax = plt.subplot(1,2,l+1)
    adjust_spines(ax, ['left', 'bottom'])
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_color('dimgrey')
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
    ax.tick_params(axis='x',labelsize=7,pad=4)
    ax.tick_params(axis='y',labelsize=7,pad=4)
    
    ax.yaxis.grid(zorder=1,color='darkgrey',alpha=0.6,clip_on=False)
    
    plt.axhline(y=0,color='dimgrey',linestyle='-',linewidth=2)
    plt.plot(years[:-length[l]+1],linetrend,linestyle='--',color='teal',label=r'\textbf{ERA5}',
             linewidth=2.5,clip_on=False,zorder=3,dashes=(1,0.3))
    plt.plot(years2[:-length[l]+1],linetrendc,linestyle='-',color='maroon',label=r'\textbf{GISTEMPv4}',
             linewidth=2.5,clip_on=False,zorder=2)
            
    if l == 0:
        plt.ylabel(r'\textbf{10-Year Trends in TAS [$^{\circ}$C/decade]}',fontsize=10,color='k')
        plt.yticks(np.arange(-10,11,1),map(str,np.round(np.arange(-10,11,1),2)),size=9)
        plt.xticks(np.arange(1900,2101,15),map(str,np.arange(1900,2101,15)),size=9)
        plt.xlim([1900,2010])   
        plt.ylim([-5,5])
        
        plt.text(1900,4.432,r'\textbf{[a]}',color='k',fontsize=8)
    elif l ==1:
        plt.ylabel(r'\textbf{30-Year Trends in TAS [$^{\circ}$C/decade]}',fontsize=10,color='k')
        plt.yticks(np.arange(-10,11,0.2),map(str,np.round(np.arange(-10,11,0.2),2)),size=9)
        plt.xticks(np.arange(1900,2101,15),map(str,np.arange(1900,2101,15)),size=9)
        plt.xlim([1900,2010])   
        plt.ylim([-1,1])
        
        plt.text(1900,0.88,r'\textbf{[b]}',color='k',fontsize=8)

plt.tight_layout()
leg = plt.legend(shadow=False,fontsize=20,loc='upper center',
            bbox_to_anchor=(-0.1, 1.17),fancybox=True,ncol=15,frameon=False,
            handlelength=1,handletextpad=0.5)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
plt.subplots_adjust(top=0.9)
plt.savefig(directoryfigure + 'GMST_TrendsComparison.png',dpi=300)