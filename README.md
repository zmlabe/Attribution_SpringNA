# Attribution_SpringNA
Climate drivers of the springtime North America cooling pattern

###### Under construction... ```[Python 3.9]```

## Contact
Zachary Labe - [Research Website](http://sites.uci.edu/zlabe/) - [@ZLabe](https://twitter.com/ZLabe)

## Description
+ ```Scripts/```: Main [Python](https://www.python.org/) scripts/functions used in data analysis and plotting
+ ```requirements.txt```: List of environments and modules associated with the most recent version of this project. A Python [Anaconda3 Distribution](https://docs.continuum.io/anaconda/) was used for our analysis. Tools including [NCL](https://www.ncl.ucar.edu/), [CDO](https://code.mpimet.mpg.de/projects/cdo), and [NCO](http://nco.sourceforge.net/) were also used for initial data processing.

## Data
+ Berkeley Earth Surface Temperature project (BEST) : [[DATA]](http://berkeleyearth.org/data/)
    + Rohde, R. and Coauthors (2013) Berkeley earth temperature averaging process. Geoinform Geostat Overv. doi:10.4172/2327-4581.1000103 [[PUBLICATION]](http://www.scitechnol.com/2327-4581/2327-4581-1-103.php)
+ CESM1 Large Ensemble Project (LENS1) : [[DATA]](http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html)
    + Kay, J. E and Coauthors, 2015: The Community Earth System Model (CESM) large ensemble project: A community resource for studying climate change in the presence of internal climate variability. Bull. Amer. Meteor. Soc., 96, 1333–1349, doi:10.1175/BAMS-D-13-00255.1 [[PUBLICATION]](http://journals.ametsoc.org/doi/full/10.1175/BAMS-D-13-00255.1)
+ CESM2 Large Ensemble Project (LENS2) : [[DATA]](https://www.cesm.ucar.edu/projects/community-projects/LENS2/)
    + Rodgers, K. B., Lee, S. S., Rosenbloom, N., Timmermann, A., Danabasoglu, G., Deser, C., ... & Yeager, S. G. (2021). Ubiquity of human-induced changes in climate variability. Earth System Dynamics Discussions, 1-22, doi:10.1175/BAMS-D-13-00255.1 [[PUBLICATION]](https://esd.copernicus.org/preprints/esd-2021-50/)
+ ERA5 : [[DATA]](https://cds.climate.copernicus.eu/cdsapp#!/home)
    + Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., ... & Simmons, A. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, doi:10.1002/qj.3803 [[PUBLICATION]](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3803)
+ GISTEMPv4 : [[DATA]](https://data.giss.nasa.gov/gistemp/)
    + Lenssen, N., G. Schmidt, J. Hansen, M. Menne, A. Persin, R. Ruedy, and D. Zyss, 2019: Improvements in the GISTEMP uncertainty model. J. Geophys. Res. Atmos., 124, no. 12, 6307-6326, doi:10.1029/2018JD029522.[[PUBLICATION]](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018JD029522)
+ GFDL CM3 Large Ensemble : [[DATA]](https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.CLIVAR_LE.html)
    + Sun, L., Alexander, M., & Deser, C. (2018). Evolution of the global coupled climate response to Arctic sea ice loss during 1990–2090 and its
contribution to climate change. Journal of Climate, 31(19), 7823–7843. doi:10.1175/JCLI-D-18-0134.1 [[PUBLICATION]](https://journals.ametsoc.org/view/journals/clim/31/19/jcli-d-18-0134.1.xml)
+ GFDL ES2M Large Ensemble : [[DATA]](https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.CLIVAR_LE.html)
    + Rodgers, K. B., Lin, J., & Frölicher, T. L. (2015). Emergence of multiple ocean ecosystem drivers in a large ensemble suite with an Earth system. doi:10.5194/bg-12-3301-2015 [[PUBLICATION]](https://bg.copernicus.org/articles/12/3301/2015/)
model. Biogeosciences, 12(11), 3301–3320. https://doi.org/10.5194/BG-12-3301-2015 
+ GFDL FLOR: Forecast-oriented Low Ocean Resolution : [[DATA]](https://www.gfdl.noaa.gov/cm2-5-and-flor/)
    + Vecchi, G. A., Delworth, T., Gudgel, R., Kapnick, S., Rosati, A., Wittenberg, A. T., ... & Zhang, S. (2014). On the seasonal forecasting of regional tropical cyclone activity. Journal of Climate, 27(21), 7994-8016. doi:10.1175/JCLI-D-14-00158.1 [[PUBLICATION]](https://journals.ametsoc.org/view/journals/clim/27/21/jcli-d-14-00158.1.xml)
+ GFDL SPEAR: Seamless System for Prediction and EArth System Research : [[DATA]](https://www.gfdl.noaa.gov/spear_large_ensembles/)
    + Delworth, T. L., Cooke, W. F., Adcroft, A., Bushuk, M., Chen, J. H., Dunne, K. A., ... & Zhao, M. (2020). SPEAR: The next generation GFDL modeling system for seasonal to multidecadal prediction and projection. Journal of Advances in Modeling Earth Systems, 12(3), e2019MS001895. doi:10.1029/2019MS001895 [[PUBLICATION]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019MS001895)
+ MPI-ESM1.2LR-LE [[DATA]](https://esgf-node.llnl.gov/search/cmip6/)
    + Mauritsen, T., Bader, J., Becker, T., Behrens, J., Bittner, M., Brokopf, R., ... & Roeckner, E. (2019). Developments in the MPI‐M Earth System Model version 1.2 (MPI‐ESM1. 2) and its response to increasing CO2. Journal of Advances in Modeling Earth Systems, 11(4), 998-1038. doi:10.1029/2018MS001400 [[PUBLICATION]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001400)
+ Multi-Model Large Ensemble (SMILE) : [[DATA]](https://www.cesm.ucar.edu/projects/community-projects/MMLEA/)
    + Deser, C., Phillips, A. S., Simpson, I. R., Rosenbloom, N., Coleman, D., Lehner, F., ... & Stevenson, S. (2020). Deser, C., Lehner, F., Rodgers, K. B., Ault, T., Delworth, T. L., DiNezio, P. N., ... & Ting, M. (2020). Insights from Earth system model initial-condition large ensembles and future prospects. Nature Climate Change, 1-10. doi:10.1038/s41558-020-0731-2 [[PUBLICATION]](https://www.nature.com/articles/s41558-020-0731-2)
+ Sixth Version of Model for Interdisciplinary Research on Climate (MIROC6) : [[DATA]](https://climexp.knmi.nl/selectfield_cmip6.cgi?id=someone@somewhere)
    + Tatebe, H., Ogura, T., Nitta, T., Komuro, Y., Ogochi, K., Takemura, T., ... & Kimoto, M. (2019). Description and basic evaluation of simulated mean state, internal variability, and climate sensitivity in MIROC6. Geoscientific Model Development, 12(7), 2727-2765. doi:10.5194/gmd-12-2727-2019 [[PUBLICATION]](https://gmd.copernicus.org/articles/12/2727/2019/gmd-12-2727-2019.html)
+ NCEP/DOE Reanalysis II : [[DATA]](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html)
    + NCEP-DOE AMIP-II Reanalysis (R-2): M. Kanamitsu, W. Ebisuzaki, J. Woollen, S-K Yang, J.J. Hnilo, M. Fiorino, and G. L. Potter. 1631-1643, Nov 2002, Bulletin of the American Meteorological Society. doi:10.1175/BAMS-83-11-1631 [[PUBLICATION]](https://journals.ametsoc.org/view/journals/bams/83/11/bams-83-11-1631.xml)
+ NOAA-CIRES-DOE Twentieth Century Reanalysis (20CRv3) : [[DATA]](https://psl.noaa.gov/data/gridded/data.20thC_ReanV3.html)
    + Slivinski, L. C., Compo, G. P., Whitaker, J. S., Sardeshmukh, P. D., Giese, B. S., McColl, C., ... & Wyszyński, P. (2019). Towards a more reliable historical reanalysis: Improvements for version 3 of the Twentieth Century Reanalysis system. Quarterly Journal of the Royal Meteorological Society, 145(724), 2876-2908. doi:10.1002/qj.3598 [[PUBLICATION]](https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3598)
+ NOAA Merged Land Ocean Global Surface Temperature Analysis (NOAAGlobalTemp v5.1) : [[DATA]](https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v5.1/access/gridded/)
    + Vose, R.S., B. Huang, X. Yin., D. Arndt, D.R. Easterling, J.H. Lawrimore, M.J. Menne, A. Sanchez-Lugo, H.-M. Zhang, 2021: Implementing full spatial coverage in NOAA’s global temperature analysis. Geophysical Research Letters, 48, e2020GL090873. doi:10.1029/2020GL090873 [[PUBLICATION]](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL090873)
+ SMHI Large Ensemble (SMHI-LENS) : [[DATA]](https://esgf-node.llnl.gov/search/cmip6/)
    + Wyser, K., Koenigk, T., Fladrich, U., Fuentes-Franco, R., Karami, M. P., & Kruschke, T. (2021). The SMHI large ensemble (SMHI-lens) with EC-Earth3. 3.1. Geoscientific Model Development, 14(7), 4781-4796. doi:10.5194/gmd-14-4781-2021 [[PUBLICATION]](https://gmd.copernicus.org/articles/14/4781/2021/)

## Publications


## Conferences/Presentations
+ **[2]** **Labe, Z.M.**, N.C. Johnson, and T.L Delworth. Climate drivers of the recent springtime cooling pattern in northern North America, 36th Conference on Climate Variability and Change, Virtual Attendance (Jan 2023). [[Abstract]](https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Paper/415409)[[Slides]](https://www.slideshare.net/ZacharyLabe/climate-drivers-of-the-recent-springtime-cooling-pattern-in-northern-north-america)
+ **[1]** **Labe, Z.M.**, N.C. Johnson, and T.L Delworth. Identifying the drivers of the observed springtime cooling trend in northern North America with large ensemble simulations, 2022 American Geophysical Union Annual Meeting, Chicago, IL (Dec 2022). [[Abstract]](https://agu.confex.com/agu/fm22/meetingapp.cgi/Paper/1111909)[[Poster]](https://zacklabe.files.wordpress.com/2022/12/labejohnsondelworth_agu_largeensembles2022_poster.pdf)
