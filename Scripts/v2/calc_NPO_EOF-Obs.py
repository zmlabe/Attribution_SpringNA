"""
Calculate Alaska Blocking Index

Reference : Ballinger et al. 2022 (IJOC)
Author    : Zachary M. Labe
Date      : 16 November 2022
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
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/v2/obs/NPO/' 
directorydata = '/work/Zachary.Labe/Data/' 

### Create months to loop through
sliceperiod = 'FMA'
yearsrange = np.arange(1979,2020+1,1)

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
timelength = 15
eofall = []
pcall = []
varianceall = []
for i,yr in enumerate(range(yearsrange.min(),yearsrange.max()-timelength+1,1)):
    print (i,yr)
    variq = 'SLP'
    yearsall = np.arange(1979,2021+1,1)
    sliceshape = 3
    slicenan = 'nan'
    addclimo = True
    yearmin = yr
    yearmax = yr+timelength

    ### Read data
    latobs,lonobs,lev,varn = ERA.read_ERA5_monthly1x1(variq,directorydata,sliceperiod,
                                                  yearsall,sliceshape,addclimo,slicenan,'surface')
    
    ### Read only 1979-2020
    yearq = np.where((yearsall >= 1979) & (yearsall <= 2020))[0]
    years = yearsall[yearq]
    var = varn[yearq,:,:]
    
    ### Calculate anomalies
    anoms = calc_anomalies(years,var)
    
    ### Detrend data
    vardt = DT.detrendDataR(anoms,'surface','monthly')
    yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
    vardt = vardt[yearq]

    ### Calculate North Pacific Box
    lon2,lat2 = np.meshgrid(lonobs,latobs)
    lonq = np.where((lonobs >=120) & (lonobs <=240))[0]
    latq = np.where((latobs >=20) & (latobs <=60))[0]
    lat2sq = lat2[latq,:]
    lat2s = lat2sq[:,lonq]
    lon2sq = lon2[latq,:]
    lon2s = lon2sq[:,lonq]
    lat1s = latobs[latq]
    lon1s = lonobs[lonq]
    
    anomlon = vardt[:,:,lonq]
    anoms = anomlon[:,latq,:]
    
    ### Calculate EOF (NPO is EOF2 of North Pacific)
    coslat = np.cos(np.deg2rad(lat1s)).clip(0.,1.)
    wgts = np.sqrt(coslat)[...,np.newaxis]
    solver = Eof(anoms,weights=wgts)
    
    ### First 2 eofs
    eofpattern = solver.eofs(neofs=2)
    eof1 = eofpattern[0]
    eof2 = eofpattern[1]
    
    ### Calculate and standardize the PCs
    pc = solver.pcs(npcs=2,pcscaling=0) # not scaled
    pc1 = sts.zscore(pc[:,0])
    pc2 = sts.zscore(pc[:,1])
    
    ### Variance explained by eof pattern
    variance_fraction = solver.varianceFraction(neigs=2)
    variance1 = variance_fraction[0]
    variance2 = variance_fraction[1]
    
    ### Make same sign
    if eof2[2,70] > 0:
        eof2 = eof2 * -1
        pc2 = pc2 * -1
    if np.nanmean(eof1) < 0:
        eof1 = eof1 * -1
        pc1 = pc1 * -1
    
    eofall.append(eof2)
    pcall.append(pc2)
    varianceall.append(variance2)
    
### Calculate location of max node for EOF2
maxlat = []
maxlon = []
for d in range(len(eofall)):
    latwhere = np.unravel_index(eofall[d].argmax(),eofall[d].shape)[0]
    lonwhere = np.unravel_index(eofall[d].argmax(),eofall[d].shape)[1]
    latnode = lat1s[latwhere]
    lonnode = lon1s[lonwhere]
    
    maxlat.append(latnode)
    maxlon.append(lonnode)
    
### Calculate location of min node for EOF2
minlat = []
minlon = []
for d in range(len(eofall)):
    latwheremin = np.unravel_index(eofall[d].argmin(),eofall[d].shape)[0]
    lonwheremin = np.unravel_index(eofall[d].argmin(),eofall[d].shape)[1]
    latnodemin = lat1s[latwheremin]
    lonnodemin = lon1s[lonwheremin]
    
    minlat.append(latnodemin)
    minlon.append(lonnodemin)
    
m = Basemap(projection='ortho',lon_0=180,lat_0=55,resolution='l',area_thresh=10000)
m.drawcoastlines(color='dimgrey',linewidth=1)
m.drawstates(color='dimgrey',linewidth=0.5)
m.drawcountries(color='dimgrey',linewidth=0.5)   
cs1 = m.contourf(lon2s,lat2s,eof2,np.arange(-0.03,0.031,0.005),extend='both',latlon=True,cmap=cmocean.cm.balance)

yearslist = ['1980-1994','1985-1999','1990-2004','1995-2009','2000-2014','2005-2019']
eofpatterns_all = [eofall[1],eofall[6],eofall[11],eofall[16],eofall[21],eofall[26]]
variance_all = [varianceall[1],varianceall[6],varianceall[11],varianceall[16],varianceall[21],varianceall[26]]

letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"]
limit = np.arange(-0.04,0.041,0.005)
barlim = np.round(np.arange(-0.04,0.041,0.04),2)
cmap = cmocean.cm.balance
label = r'\textbf{EOF2 - NPO}'

fig = plt.figure(figsize=(5,10))
for i in range(len(eofpatterns_all)):
    plt.subplot(6,1,i+1)

    m = Basemap(projection='cyl', llcrnrlon=120, llcrnrlat=20,
                urcrnrlon=240, urcrnrlat=60,resolution='l',area_thresh=1000)
    m.drawcoastlines(color='darkgrey',linewidth=1.4)
    m.drawstates(color='darkgrey',linewidth=0.5)
    m.drawcountries(color='darkgrey',linewidth=0.5)
    
    x, y = np.meshgrid(lon1s,lat1s)
       
    circle = m.drawmapboundary(fill_color='dimgrey',color='dimgray',
                      linewidth=0.7)
    circle.set_clip_on(False)
    
    cs1 = m.contourf(x,y,eofpatterns_all[i],limit,extend='both',latlon=True)
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
                
plt.savefig(directoryfigure + 'EOF2_NPO_%s_AllTimePeriods.png' % sliceperiod,dpi=300)

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

plt.plot(years[:len(minlon)],minlon,linestyle='-',color='maroon',linewidth=2,clip_on=False,label='Nothern Node')
plt.plot(years[:len(maxlon)],maxlon,linestyle='--',dashes=(1,0.3),color='teal',linewidth=2,clip_on=False,label='Southern Node')

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
plt.title(r'\textbf{Shift in the NPO centers of action - %s}' % sliceperiod,fontsize=12,color='k')
plt.savefig(directoryfigure + 'EOF2_NPOshift_obs_%s.png' % sliceperiod,dpi=300)
