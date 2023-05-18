"""
Plot trends in FMA from 1979 to 2022 for GEOP
Author    : Zachary M. Labe
Date      : 6 June 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import calc_Utilities as UT
import sys
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ERA5/monthly/' 

### Parameters
monthq = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
monthstext = ['OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
variq = 'GEOP'
sliceperiod = 'none'
years = np.arange(1979,2021+1,1)
sliceshape = 4
slicenan = 'nan'
addclimo = True
newon = False

### Read data
data = Dataset(directorydata + 'GEOP_1979-2021.nc')
geop = data.variables['GEOP'][:]
level = data.variables['level'][:]
lat = data.variables['latitude'][:]
lon = data.variables['longitude'][:]
data.close()

### Calculate polar cap
latq = np.where((lat >= 65))[0]
latcap = lat[latq]
loncap2,latcap2 = np.meshgrid(lon,latcap)
cap = UT.calc_weightedAve(geop[:,:,latq,:],latcap2)

### Reshape for year,month,level
caps = cap.reshape(cap.shape[0]//12,12,level.shape[0])

### Calculate trends
trend = np.empty((caps.shape[1],caps.shape[2]))
for mo in range(caps.shape[1]):
    for l in range(caps.shape[2]):
        mask = np.isfinite(caps[:,mo,l])
        x = np.arange(years.shape[0])
        y = caps[:,mo,l]
        
        if np.sum(mask) == y.shape[0]:
            xx = x
            yy = y
        else:
            xx = x[mask]
            yy = y[mask]
        
        if np.isfinite(np.nanmean(yy)):
            trend[mo,l],intercepts,r_value,p_value,std_err = sts.linregress(xx[:-2],yy[:-2])
        else:
            trend[mo,l] = np.nan
            
### Change time for plotting (October - May) for decade trends
trendwinter = np.append(trend[9:,:],trend[:5,:],axis=0) * 10

###########################################################################
###########################################################################
###########################################################################
##### Plot profiles
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
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
        
### Set limits for contours and colorbars6
limit = np.arange(-50,50.1,0.1)
barlim = np.arange(-50,51,25)
cmap = cmocean.cm.balance
label = r'\textbf{ZCAP [m/decade]}'
zscale = np.array([1000,700,500,300,200,100,50,30,10])
timeq,levq = np.meshgrid(np.arange(len(monthstext)),level)
        
fig = plt.figure(figsize=(10,4))

### Create plot
ax1 = plt.subplot(111)
ax1.spines['top'].set_color('dimgrey')
ax1.spines['right'].set_color('dimgrey')
ax1.spines['bottom'].set_color('dimgrey')
ax1.spines['left'].set_color('dimgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.tick_params(axis='x',direction='out',which='major',pad=3,
            width=2,color='dimgrey')   
ax1.tick_params(axis='y',direction='out',which='major',pad=3,
            width=2,color='dimgrey')
plt.gca().axes.get_xaxis().set_visible(True)
plt.gca().axes.get_yaxis().set_visible(True)
plt.ylabel(r'\textbf{Pressure [hPa]}',color='k',fontsize=7)

ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

### Plot contours
cs = plt.contourf(timeq,levq,trendwinter.transpose(),limit,extend='both')

cs.set_cmap(cmap)

plt.gca().invert_yaxis()
plt.yscale('log')

plt.xticks(np.arange(0,7+1,1),monthstext,fontsize=4)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)

plt.xlim([0,7])
plt.ylim([1000,10])
plt.minorticks_off()

###########################################################################
plt.tight_layout()
cbar_ax = fig.add_axes([0.33,0.08,0.4,0.03])                    
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(label,fontsize=11,color='dimgrey',labelpad=1.4)  

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001,labelsize=7)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(bottom=0.17,hspace=0.08,wspace=0.08)    
plt.savefig(directoryfigure + 'Obs_Trends-ZCap_Oct-Apr_1979-2019.png',dpi=300)
print('Completed: Script done!')
