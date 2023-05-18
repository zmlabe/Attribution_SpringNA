"""
Plot EOF for NPO using SPEAR data

Author    : Zachary M. Labe
Date      : 5 December 2022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import calc_Utilities as UT
import sys
import read_ERA5_monthly1x1 as ERA
import calc_DetrendData as DT
from eofs.standard import Eof
import scipy.stats as sts
import cmocean

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/SPEAR/NPO/' 
directorydata = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/v2/SPEAR/NPO/' 

### Create months to loop through
sliceperiod = 'April'
years = np.arange(1979,2020+1,1)

### Read in data
eof = np.load(directorydata + 'EOFinfo_NPO_coupledSPEAR_%s.npz' % sliceperiod)
eofalln = eof['eofall']
variancealln = eof['varianceall']
timelength = eof['timelength']
yearsall = eof['years']

grid = np.load(directorydata + 'LatLon_NPO_coupledSPEAR_%s.npz' % sliceperiod)
latmin = grid['minall_lat']
latmax = grid['maxall_lat']
lonmin = grid['minall_lon']
lonmax = grid['maxall_lon']

matrix = np.load(directorydata + 'LatLonGrid_NPO_coupledSPEAR_%s.npz' % sliceperiod)
lat2s = matrix['lat2s']
lon2s = matrix['lon2s']

### Calculate ensemble means
minlonm = np.nanmean(lonmin,axis=0)
maxlonm = np.nanmean(lonmax,axis=0)
varianceall = np.nanmean(variancealln,axis=1)

### Ensemble mean of EOF
eofall = np.nanmean(eofalln,axis=1)*-1
meaneof = np.nanmean(eofall,axis=0)*-1
    
m = Basemap(projection='ortho',lon_0=180,lat_0=55,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)   
cs1 = m.contourf(lon2s,lat2s,meaneof,np.arange(-0.03,0.031,0.005),extend='both',latlon=True,cmap=cmocean.cm.balance)

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
eofpatterns_all = [eofall[1],eofall[6],eofall[11],eofall[16],eofall[21],eofall[26]]
variance_all = [varianceall[1],varianceall[6],varianceall[11],varianceall[16],varianceall[21],varianceall[26]]

letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-0.01,0.011,0.001)
barlim = np.round(np.arange(-0.01,0.011,0.01),2)
cmap = cmocean.cm.balance
label = r'\textbf{EOF2 - NPO -- SPEAR}'

fig = plt.figure(figsize=(5,10))
for i in range(len(eofpatterns_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=120, llcrnrlat=20,
                urcrnrlon=240, urcrnrlat=60,resolution='l',area_thresh=1000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    # x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(lon2s,lat2s,eofpatterns_all[i],limit,extend='both',latlon=True)
    cs1.set_cmap(cmap) 
    
    if i == 0:
        parallels = np.arange(-90,91,20)
        meridians = np.arange(-180,180,20)
        par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                        color='w',fontsize=4,zorder=40)
        mer=m.drawmeridians(meridians,labels=[False,False,True,False],linewidth=0.5,
                            fontsize=4,color='w',zorder=40)
    elif any([i==1,i==2,i==3,i==4]):
        parallels = np.arange(-90,91,20)
        meridians = np.arange(-180,180,20)
        par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                        color='w',fontsize=4,zorder=40)
        mer=m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.5,
                            fontsize=4,color='w',zorder=40)
    elif i == 5:
        parallels = np.arange(-90,91,20)
        meridians = np.arange(-180,180,20)
        par=m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                        color='w',fontsize=4,zorder=40)
        mer=m.drawmeridians(meridians,labels=[False,False,False,True],linewidth=0.5,
                            fontsize=4,color='w',zorder=40)
     
    plt.text(104,26,r'\textbf{%s}' % yearslist[i],rotation=90,color='k',fontsize=14)   
    plt.text(121,55,r'\textbf{[%s - %s\%%]}' % (letters[i],np.round(variance_all[i]*100,1)),fontsize=8,color='k',zorder=100)
    
###############################################################################
cbar_ax1 = fig.add_axes([0.37,0.06,0.3,0.02])                
cbar1 = fig.colorbar(cs1,cax=cbar_ax1,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar1.set_label(label,fontsize=10,color='dimgrey',labelpad=1.8)  
cbar1.set_ticks(barlim)
cbar1.set_ticklabels(list(map(str,barlim)))
cbar1.ax.tick_params(axis='x', size=.01,labelsize=7)
cbar1.outline.set_edgecolor('dimgrey')

# plt.tight_layout()
plt.subplots_adjust(hspace=-0.01)
                
plt.savefig(directoryfigure + 'EOF2_NPO_%s_AllTimePeriods_SPEAR.png' % sliceperiod,dpi=300)

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
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',
               labelbottom='off',bottom='off')
ax.yaxis.grid(zorder=2,color='darkgrey',alpha=1,clip_on=False,linewidth=0.5)

# for ens in range(len(minlonm)):
#     plt.plot(years[:len(lonmin[ens])],lonmin[ens],linestyle='-',color='teal',linewidth=0.4,clip_on=False,alpha=0.3)

plt.plot(years[:len(maxlonm)],maxlonm,linestyle='-',color='maroon',linewidth=2,clip_on=False,label='Northern Node')
plt.plot(years[:len(minlonm)],minlonm,linestyle='--',dashes=(1,0.3),color='teal',linewidth=2,clip_on=False,label='Southern Node')

leg = plt.legend(shadow=False,fontsize=11,loc='upper center',
            fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,0.1),
            labelspacing=1,columnspacing=1,handletextpad=0.4)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())

plt.xticks(np.arange(1979,2021,5),map(str,np.arange(1979,2021,5)),fontsize=7)
plt.yticks(np.arange(100,300,10),map(str,np.arange(100,300,10)),fontsize=7)
plt.xlim([1979,2020])
plt.ylim([160,240])

plt.ylabel(r'\textbf{Lngitude [$^{\circ}$]}',fontsize=8,color='k')
plt.title(r'\textbf{Shift in the NPO centers of action - %s -- SPEAR}' % sliceperiod,fontsize=12,color='k')
plt.savefig(directoryfigure + 'EOF2_NPOshift_SPEAR_%s.png' % sliceperiod,dpi=300)
