% Run all the necessari codes to make the figures of the article from the
% time series data

%% 1- Compute the wavelets of the Ciutadella sea level and atmospheric pressure time series
C1_compute_wavelets
% the data is stored in: ./data/ciutadella_spectral_data.mat
disp('C1 done')
disp('--------------------------------------------------------------------')
close all; clear all;
%% 2- Find the meteotsunami events of the period of study
% the criteria to define a meteotsunami event is explained in the section
% 2.1. of the paper
C2_select_events
disp('C2 done')
disp('--------------------------------------------------------------------')
close all; clear all;

%% 3- Cut the atmospheric pressure time series in the duration of each event and store it un a new folder
C3_cut_time_series
disp('C3 done')
disp('--------------------------------------------------------------------')
close all; clear all;

%% 4- Estimate the propagation velocities in each event
% 1 first for all the 2 - 120 minutes band
C4_prop_vel_allfreq
disp('C4a done')
% 2 nd for different frequency bands
C4_prop_vel_freq
disp('C4b done')
disp('--------------------------------------------------------------------')
close all; clear all;

%% 5- Compute the maximum lagged correlation matrix during each event
C5_lag_corr_events_12h
disp('C5 done')
disp('--------------------------------------------------------------------')
%--------------------------------------------------------------------------
close all; clear all;

%% Figures:
cd ./figure_codes/
addpath '../';
% Figure 2:
Fig2_allin1
Fig2_st_wl
close all; clear all;

% Figure 3: 
Fig3_spectra
close all; clear all;

% Figure 4:
Fig4_scatter_patm_sl
close all; clear all;

%% Figure 5 & 6:
Fig5i6_EOFs
close all; clear all;

% Figure 7: 
Fig7_vel_ST
close all; clear all;

% Figure 8: 
Fig8_vel_all
close all; clear all;



