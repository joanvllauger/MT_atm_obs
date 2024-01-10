# MT_atm_obs
This folder contains all the codes necessary to reproduce the figure from the article:

Observational characterization of atmospheric disturbances generating meteotsunamis in the Balearic Islands

Joan Villalonga*(1,2), Sebastià Monserrat (1), Damià Gomis (1,3), Gabriel Jordà*(2)

(1) Departament de Física (UIB), Palma, Spain.

(2) Centre Oceanogràfic de Balears, CN-Instituto Español de Oceanografía (IEO-CSIC), Palma, Spain. 

(3) Institut Mediterrani d’Estudis Avançats (UIB-CSIC), Esporles, Spain.

Corresponding email: joan.villalonga@uib.cat 

1- Data 

The raw data files used are too heavy to be uploaded to GitHub reason why they must be found in the Zenodo repository: 10.5281/zenodo.10478602

There are 7 data files:

Atm_pres_all: containing the atmospheric pressure time series measured in the different meteorological stations used in the work. Each station contain its name and position in coordinates. All the time series have a temporal resolution of 1 min. The data have been obtained from BalearsMeteo (http://balearsmeteo.com/) and from SOCIB (https://www.socib.es/)
ciutadella_SL_AtmPres: containing the sea level and atmospheric pressure records in Ciutadella from 2018 to 2021. All the time series have a temporal resolution of 1 min. The data have been provided by PortIB (https://www.portsib.es/ca/paginas/inici).
ciutadella_SL_long: containing the sea level records in Ciutadella from 2014 to 2021. All the time series have a temporal resolution of 1 min. The data have been provided by PortIB (https://www.portsib.es/ca/paginas/inici).
ciutadella_spectral_data: containing the sea level and atmospheric pressure power wavelet spectra in Ciutadella from 2018 to 2021. Computed from the data in ciutadella_SL_AtmPres.
corr_rissagues_1min_allfreq_12h: containing the maximum lagged correlation matrices between the atmospheric pressure time series measured at the 12h surrounding each meteotsunami event in 2021. They have been computed from the data in Atm_pres_all.
sepic_index_vars: containing the five ERA5 1-hour time series of the variables used to compute the meteotsunami index as described in Sepic, et al,. 2016.
wind_ciutadella: containing the time series of the wind speed and direction provided by ERA5 reanalysis over Ciutadella during the period of study.

Once the data have been downloaded it must be placed inside a directory called data in the main directory.

2- Before running

Before running the programs developed for this study two Matlab packages must be downloaded and stored in the main directory:

- Wavelet package:
Torrence, C. and G. P. Compo, 1998: A Practical Guide to Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
Can be downloaded from: https://paos.colorado.edu/research/wavelets/

- M_map package:
Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.
Can be downloaded from: https://www.eoas.ubc.ca/~rich/map.html

Moreover, all these scripts have been run in Matlab R2021a, the may present some problemes in other versions.

3- Data processing 

The data processing codes used before plotting the figures are placed in the main directory. 

There are 6 different codes that must be run in order:

-  C1_compute_wavelets: compute the wavelet transform for the sea level and the atmospheric pressure time series
-  C2_select_events: use the algorithm described in the methods section of the article to generate a list of meteotsunami events.
-  C3_cut_time_series: uses the meteotsunami list generated in C2 to cut the atmospheric pressure time series during each meteotsunami event and store this data in the event folder located in ./meteotsunamis/meteotsunami_yyyymmdd being yyyymmdd the date of the event.
-  C4_prop_vel_allfreq: uses the atmospheric pressure time series cut for each event in C3 to estimate the propagation velocity of the atmospheric perturbation generation in each of the events. To understand how the propagation velocity is estimated read Supplementary Information 1. This code does the propagation velocity estimation for the frequency band ranging from 2 minutes to 2 hours of period.
-  C4_prop_vel_freq: uses the atmospheric pressure time series cut for each event in C3 to estimate the propagation velocity of the atmospheric perturbation generation in each of the events. To understand how the propagation velocity is estimated read Supplementary Information 1. This code does the propagation velocity estimation for the different frequency bands that can be found in Supplementary Information 5.
-  C5_lag_corr_events_12h: This code generates the maximum lagged correlation coefficient matrix between the atmospheric pressure time series measured at all the meteorological stations during the 12 hours centered at the time of each meteotsunami event. These correlation matrices will be stored and later used to estimate the EOFs of the atmospheric disturbances.

** Note that by combining this repository with the data in the Zenodo XXXXX the reader have already the output of running all these functions. Therefore it should not be necessary to run them before using the scripts that generate the figure. 

4- Figures

Once all the 6 programs have been successfully run the reader can go to the figure_codes to generate any of the figures that are in the manuscript. 
