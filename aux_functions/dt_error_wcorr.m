function [err]=dt_error_wcorr(x)
global corr_m gamma_i Dt di sig_c
%--- Parametres
c=x(1);
gamma=x(2);
%--- Importam les dades observades (Dt, gamma_i,di)
%load observed_data.mat
%--- Calculem el pes
w=exp(-(1-corr_m).^2/(2*sig_c^2));
%--- Calculem el dt teòric
dt=di/c.*cos(gamma_i-gamma);
%--- Calculem l'error
err=sum(w.*(Dt-dt).^2);
end

