"""
Plots Arctic sea ice extent decadal trends bar graph (sea ice index)

Website   : https://nsidc.org/data/seaice_index
Author    : Zachary M. Labe
Date      : 8 October 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import datetime

### Directory and time
directoryfigure = '/Users/zlabe/Documents/SciComm/SIE_DecadalTrends/'
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr

trends = np.array([-3.10,-2.87,-2.58,-2.65,-2.66,
                   -3.99,-7.45,-10.42,-12.66,
                   -9.84,-4.97,-3.50])

### Make plot
matplotlib.rc('savefig', facecolor='black')
matplotlib.rc('axes', edgecolor='darkgrey')
matplotlib.rc('xtick', color='darkgrey')
matplotlib.rc('ytick', color='darkgrey')
matplotlib.rc('axes', labelcolor='darkgrey')
matplotlib.rc('axes', facecolor='black')
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})

fig = plt.figure()
ax = plt.subplot(111) 

### Adjust axes in time series plots 
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

ax.yaxis.grid(zorder=1,color='w',alpha=0.15)

recdiff_masked = np.ma.masked_less_equal(trends, 0)

rects1 = plt.bar(np.arange(len(trends))+0.5,trends,color='crimson',
        edgecolor='crimson',zorder=9) 
rects2 = plt.bar(np.arange(len(trends))+0.5,recdiff_masked,
        color='deepskyblue',edgecolor='deepskyblue',zorder=9) 

### Define date
xlabels = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',
          r'Aug',r'Sep',r'Oct',r'Nov',r'Dec',r'Jan'] 

plt.text(0,-12.6,r'\textbf{DATA:} National Snow \& Ice Data Center, Boulder CO [1979-2020; Sea Ice Index v3]',
         fontsize=5,rotation='horizontal',ha='left',color='darkgrey')
plt.text(0,-13.2,r'\textbf{SOURCE:} ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/',
         fontsize=5,rotation='horizontal',ha='left',color='darkgrey')
plt.text(0,-13.8,r'\textbf{GRAPHIC:} Zachary Labe (@ZLabe)',
         fontsize=5,rotation='horizontal',ha='left',color='darkgrey')   
            
           
adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(0)
ax.spines['left'].set_linewidth(2)
ax.tick_params(axis='both', direction='out',length=5.5,width=2,
               which='major',pad=3,labelcolor='darkgrey')

plt.ylabel(r'\textbf{\% / decade}',
           fontsize=8,alpha=1,color='w',labelpad=0.5) 

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2.,1,
                r'\textbf{%2.2f}' % height,ha='center',va='top',zorder=11,color='w',
                fontsize=7)
autolabel(rects1)

plt.xticks(np.arange(0,12,1)+0.5,xlabels,rotation=0,fontsize=6)
ylabels = map(str,np.arange(-14,16,2))
plt.yticks(np.arange(-14,16,2),ylabels,fontsize=6)
plt.ylim([-14,14])
plt.xlim([0,12])

plt.title(r'\textbf{ARCTIC SEA ICE EXTENT TRENDS}',fontsize=23,color='w',alpha=1) 

plt.savefig(directoryfigure + 'NSIDC_SIE_DecadalTrends.png',dpi=900)
