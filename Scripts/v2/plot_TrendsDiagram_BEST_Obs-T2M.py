"""
Plot trend snapshot from Berkeley Earth

Author    : Zachary M. Labe
Date      : 4 January 2023
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cmocean
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import calc_Utilities as UT
import sys
from netCDF4 import Dataset
import scipy.stats as sts
import sys
import palettable.cubehelix as cm

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/Presentations/' 

### New functions
def calc_regionalAve(regionn,lat1,lon1,data):
    if regionn == 'CANMID':
        la1 = 43
        la2 = 60
        lo1 = 240
        lo2 = 295
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
    elif regionn == 'ARCTIC':
        la1 = 65
        la2 = 90
        lo1 = 0
        lo2 = 360
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
    elif regionn == 'GLOBAL':
        la1 = -90
        la2 = 90
        lo1 = 0
        lo2 = 360
        lat1q = np.where((lat1 >= la1) & (lat1 <= la2))[0]
        lon1q = np.where((lon1 >= lo1) & (lon1 <= lo2))[0]
        
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
    mean = UT.calc_weightedAve(meanbox,lat2q)

    return mean,lat1a,lat2q,lon1a,lon2q

### Parameters
monthq = ['JAN','FEB','MAR','ARP','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
datalabels = ['ERA5','BEST']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'T2M'
sliceperiod = 'April'
yearsb = np.arange(1850,2021+1,1)
years = np.arange(1850,2020+1,1)
yearse = np.arange(1950,2021+1,1)
yearsee = np.arange(1950,2020+1,1)
sliceshape = 3
slicenan = 'nan'
addclimo = True

### Read data
directorydata = directory = '/work/Zachary.Labe/Data/BEST/'
filename = 'T2M_BEST_1850-2021.nc'
data = Dataset(directory + filename,'r')
lat1 = data.variables['latitude'][:]
lon1 = data.variables['longitude'][:]
anom = data.variables['T2M'][:,:,:]
data.close()

lon2,lat2 = np.meshgrid(lon1,lat1)

### Calculate month
alldata = anom.reshape(yearsb.shape[0],12,lat1.shape[0],lon1.shape[0])
month21 = np.nanmean(alldata[:,3:4,:,:],axis=1)
month = month21[:-1,:,:]

### Calculate region
obs_reg,latr1,latr2,lonr1,lonr2 = calc_regionalAve('CANMID',lat1,lon1,month)

### Calculate trends
trend = np.full((len(obs_reg)),np.nan)
stat = np.full((len(obs_reg)),np.nan)
for yrall in range(len(obs_reg)-9):
    x = np.arange(0,len(obs_reg)-yrall)
    y = obs_reg[yrall:len(obs_reg)]

    slope, intercept, r, p, se = sts.linregress(x,y)
    
    trend[yrall] = slope * 10.
    stat[yrall] = p
    
### Calculate region
obs_regg,latr1g,latr2g,lonr1g,lonr2g = calc_regionalAve('GLOBAL',lat1,lon1,month)

### Calculate trends
trendg = np.full((len(obs_regg)),np.nan)
statg = np.full((len(obs_regg)),np.nan)
for yrall in range(len(obs_regg)-9):
    xg = np.arange(0,len(obs_regg)-yrall)
    yg = obs_regg[yrall:len(obs_regg)]

    slopeg, interceptg, rg, pg, seg = sts.linregress(xg,yg)
    
    trendg[yrall] = slopeg * 10.
    statg[yrall] = pg
    
### Calculate region
obs_rega,latr1a,latr2a,lonr1a,lonr2a = calc_regionalAve('ARCTIC',lat1,lon1,month)

### Calculate trends
trenda = np.full((len(obs_rega)),np.nan)
stata = np.full((len(obs_rega)),np.nan)
for yrall in range(len(obs_rega)-9):
    xa = np.arange(0,len(obs_rega)-yrall)
    ya = obs_rega[yrall:len(obs_rega)]

    slopea, intercepta, ra, pa, sea = sts.linregress(xa,ya)
    
    trenda[yrall] = slopea * 10.
    stata[yrall] = pa
    
# ### Read data
# directorydatae = '/work/Zachary.Labe/Data/ERA5_1x1/'
# filenamee = 'T2M_1950-2021.nc'
# datae = Dataset(directorydatae + filenamee,'r')
# lat1e = datae.variables['latitude'][:]
# lon1e = datae.variables['longitude'][:]
# anome = datae.variables['T2M'][:,:,:] - 273.15
# datae.close()

# lon2e,lat2e = np.meshgrid(lon1e,lat1e)

# ### Calculate month
# alldatae = anome.reshape(yearse.shape[0],12,lat1e.shape[0],lon1e.shape[0])
# month21e = np.nanmean(alldatae[:,3:4,:,:],axis=1)
# monthe = month21e[:-1,:,:]

# ### Calculate region
# obs_rege,latr1e,latr2e,lonr1e,lonr2e = calc_regionalAve('CANMID',lat1e,lon1e,monthe)

# ### Calculate trends
# trende = np.full((len(obs_rege)),np.nan)
# state = np.full((len(obs_rege)),np.nan)
# for yrall in range(len(obs_rege)-9):
#     xe = np.arange(0,len(obs_rege)-yrall)
#     ye = obs_reg[yrall:len(obs_rege)]
    
#     slopee, intercepte, re, pe, see = sts.linregress(xe,ye)
    
#     trende[yrall] = slopee * 10.
#     state[yrall] = p
    
# ### Combine for empty array
# empty = np.full(len(trend)-len(trende),np.nan)
# trendee = np.append(empty,trende,axis=0)

# trend = np.full((len(obs_reg),len(obs_reg)),np.nan)
# stat = np.full((len(obs_reg),len(obs_reg)),np.nan)
# for i in range(len(obs_reg)-9):
#     for j in range(9,len(obs_reg)-9):
#         trendlength = len(obs_reg) - i - j
        
#         if trendlength >= 10:
#             x = np.arange(trendlength)
#             y = obs_reg[i:trendlength+i]
            
#             slope, intercept, r, p, se = sts.linregress(x,y)
#             trend[trendlength,j] = slope * 10.
#             stat[trendlength,j] = p
#         else:
#             trend[trendlength,j] = np.nan
#             stat[trendlength,j] = np.nan
# plt.contourf(trend,np.arange(-1,1.01,0.01),extend='both',cmap=cmocean.cm.balance)

###############################################################################
###############################################################################
###############################################################################
### Plot monthly indices
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

fig = plt.figure()
ax = plt.subplot(111)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey',
               labelbottom='off',bottom='off')

plt.axhline(y=0,color='darkgrey',linewidth=2)
plt.plot(years,trend,linestyle='--',linewidth=4,color='deepskyblue',
         label=r'Region - North America',dashes=(1,0.4))
plt.plot(years,trenda,linestyle='-',linewidth=3,color='crimson',
         label=r'Arctic')
plt.plot(years,trendg,linestyle='-',linewidth=2,color='w',
         label=r'Global')

leg = plt.legend(shadow=False,fontsize=10,loc='upper center',
            fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.43,0.1),
            labelspacing=1,columnspacing=1,handletextpad=0.5)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())


plt.xticks(np.arange(1850,2021,10),map(str,np.flip(np.arange(0,172,10))),fontsize=7)
plt.yticks(np.arange(-4,2.1,0.25),map(str,np.round(np.arange(-4,2.1,0.25),2)),fontsize=7)
plt.xlim([1850,2020])
plt.ylim([-2.5,2])
plt.xlabel(r'\textbf{Length of Trend [years]}',fontsize=11,color='w')
plt.ylabel(r'\textbf{T2M [$^{\circ}$C/decade]}',fontsize=11,color='w')

plt.title(r'\textbf{TREND PERIODS IN APRIL SINCE 1850}',fontsize=20,color='w')

plt.tight_layout()
plt.savefig(directoryfigure + 'PRESENT_TrendSegmentsObs.png',dpi=300)
