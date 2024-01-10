This folder conatins all the codes necessaries to reproduce the figure from the article:

Observational characterization of atmospheric disturbances generating meteotsunamis in the Balearic Islands
Joan Villalonga*(1,2), Sebastià Monserrat (1), Damià Gomis (1,3), Gabriel Jordà*(2)
(1) Departament de Física (UIB), Palma, Spain.
(2) Centre Oceanogràfic de Balears, CN-Instituto Español de Oceanografía (IEO-CSIC), Palma, Spain. 
(3) Institut Mediterrani d’Estudis Avançats (UIB-CSIC), Esporles, Spain.
Correponding Email: joan.villalonga@uib.cat 

1- Data 
The raw data files used are too heavy to be uploaded to gighub reason why the must be found in the Zenodo repository: 

Once the data have been downloaded it must be plased inside a directory called data in the main directory.

2- Before runing
Before reuning the programs developed for this study two matlab packages must be downloaded and stored in the main direcory:

- Wavelet package:
Torrence, C. and G. P. Compo, 1998: A Practical Guide to Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
Can be downloaded from: https://paos.colorado.edu/research/wavelets/

- M_map package:
Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.
Can be downloaded from: https://www.eoas.ubc.ca/~rich/map.html

Moreover, all these scripts have been run in Matlab R2021a, the may present some problemes in other versions. 

3- Data processing 
The data processing codes used before ploting the figures are placed in the main directory. 

There are the 6 diferent codes that must be runned in order:

-  C1_compute_wavelets: compute the wavelet tranform for the sea level and the atmospheric pressure time series
-  C2_select_events: use the algorithm described in the mathods section of the article to generate a list of meteotsunami events.
-  C3_cut_time_series: uses the meteotsunami list generated in C2 to cut the atmospheric pressure time series during each meteotsunami event and store this data in the event folder located in: ./meteotsunamis/meteotsunami_yyyymmdd being yyyymmdd the date of the event.
-  C4_prop_vel_allfreq: uses the atmospheric pressure time series cut for each event in C3 to estimate the propagation velocity of the atmospheric perturbation generation each of the events. To undestand how the propagation velocity is estimated read Supplementary Information 1. This code does the propagation velocity estimation for the frequency band ranging from 2 minutes to 2 hours of period.
-  C4_prop_vel_freq: uses the atmospheric pressure time series cut for each event in C3 to estimate the propagation velocity of the atmospheric perturbation generation each of the events. To undestand how the propagation velocity is estimated read Supplementary Information 1. This code does the propagation velocity estimation for the different frequency bands that can be found in Supplementary Information 5.
-  C5_lag_corr_events_12h: This code generates the maximum lagged correlation coefficient matrix between the atmospheric pressure time series measured at all the meteotsological stations during the 12 hours centered at the time of each meteotsunami event. This correlation matrices will be stored and later used to estimate the EOFs of the atmospheric disturbances.

** Note that by combining these repository with the data in the Zenodo XXXXX the reader have already the output of running all this functions. Therefore it should not be necessary to run them before using the scripts that generate the figure. 

4- Figures
Once all the 6 programs have been succesfully run the reader can go to the figure_codes to generate any of the figures that are in the manuscript. 

For any doubts fell free to contact the authors at: joan.villalonga@uib.cat

Please acknowledge the use of the above software in any publications




