"""
Plot monthly trends in the Arctic

Author    : Zachary M. Labe
Date      : 27 June 022
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import calc_Utilities as UT
import sys
import scipy.stats as sts

### Read in data files from server
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
directoryfigure = '/home/Zachary.Labe/Research/Attribution_SpringNA/Figures/' 
directorydata = '/work/Zachary.Labe/Data/ClimateIndices/' 
directorydata2 = '/home/Zachary.Labe/Research/Attribution_SpringNA/Data/' 

monq = 'April'
yearsall = np.arange(1950,2021+1,1)
yearmin = 1979
yearmax = 2019
yearq = np.where((yearsall >= yearmin) & (yearsall <= yearmax))[0]
years = yearsall[yearq]

### Read in observations
obs = np.genfromtxt(directorydata2 + 'ERA5_T2M_BoxOfInterest_TimeSeries-%s.txt' % monq,unpack=True)
obsz = sts.zscore(obs)

### Read in teleconnections
naod = np.genfromtxt(directorydata + 'NAO/NAO_CPC_1950-2021.txt',unpack=True)
years_nao = naod[0,:]
nao_index = naod[1:,yearq]

pnad = np.genfromtxt(directorydata + 'PNA/PNA_CPC_1950-2021.txt',unpack=True)
years_pna = pnad[0,:]
pna_index = pnad[1:,yearq]

### Read in new indices
PNAd = np.genfromtxt(directorydata + 'PNA/PNAmodified_%s_%s_%s-%s_detrended.txt' % ('Z500',monq,yearmin,yearmax),unpack=True)
years_PNA = PNAd[0,:]
PNAm = sts.zscore(PNAd[1,:])

NAOd = np.genfromtxt(directorydata + 'NAO/NAOmodified_%s_%s_%s-%s_detrended.txt' % ('SLP',monq,yearmin,yearmax),unpack=True)
years_NAO = NAOd[0,:]
NAOm = sts.zscore(NAOd[1,:])

### Slice months
if monq == 'FMA':
    nao = np.nanmean(nao_index[1:4,:],axis=0)
    pna = np.nanmean(pna_index[1:4,:],axis=0)
elif monq == 'April':
    pna = np.nanmean(pna_index[3:4,:],axis=0)
    nao = np.nanmean(nao_index[3:4,:],axis=0)
elif monq == 'none':
    pna = pna_index
    nao = nao_index
else:
    print(ValueError('wrong months selected!'))
    sys.exit()
    
### Calculate trends
slope_nao,intercept_nao,r_nao,p_nao,se_nao = sts.linregress(years,nao)
slope_pna,intercept_pna,r_pna,p_pna,se_pna = sts.linregress(years,pna)

### Calculate correlations
corr_nao,pcorr_nao = sts.pearsonr(obsz,nao)
corr_pna,pcorr_pna = sts.pearsonr(obsz,pna)
