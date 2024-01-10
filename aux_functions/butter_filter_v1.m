% Filtre de IIR (butter) zero-phase filter
% Joan Villalonga 15/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%       x: original signal
%       dt: period sampling
%       T_cut: period of cut
%       order: Order of the filter
%       FiltType: can be 'low' or 'high' (in frequency)
%                 'high'--> low periods pass
%                 'low' --> 
% OUTPUS
%       y: filterred signal
%
function [y]=butter_filter_v1(x,dt,T_cut,order,FiltType)
%--- Calculem la freuqnecia de tall
Wn=(1/T_cut)*dt*2;

%--- Coefficients del filtre
[b,a]=butter(order,Wn,FiltType);
y=filtfilt(b,a,x);
end




