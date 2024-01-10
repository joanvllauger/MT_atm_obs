addpath './wavelet_package/'
addpath './aux_functions/'
% compute the wavelet transform
load('./data/ciutadella_SL_AtmPres.mat','time','SL','Patm','lon','lat','mareografs')
SL=QC_SL_v1(SL);
Patm=QC_SL_v1(Patm);
nans=isnan(SL)|isnan(Patm);

t_max_gap=time(end);
SL=rellenar_huecos(time,SL,t_max_gap,'lineal');
Patm=rellenar_huecos(time,Patm,t_max_gap,'lineal');
[wsl,period,coi,sig_sl]=wavelet_v1(SL,mode(diff(time))*24*60);
[wpa,period,coi,sig_pa]=wavelet_v1(Patm,mode(diff(time))*24*60);

save('./data/ciutadella_spectral_data.mat','wpa','period','sig_pa','wsl','coi','sig_sl','SL','time','Patm','nans')