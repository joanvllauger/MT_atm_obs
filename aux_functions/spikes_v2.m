function [y]=spikes_v2(x,k_window,n_std)
% Control de Qualitat
y=x;
% eliminam spikes
mv=movmean(y,k_window,'omitnan');
aux=abs(y-mv);
stdaux=movstd(y,k_window,'omitnan');
y(aux>stdaux*n_std)=NaN;