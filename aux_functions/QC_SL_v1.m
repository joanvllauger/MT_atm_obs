function SL=QC_SL_v1(SL)

%  figure
%  plot(SL)
%  hold on
% més enllà de 8 sigmas de tota la serie
n_sd=8;
sd=std(SL,'omitnan');
aux=abs(SL-mean(SL,'omitnan'));
SL(aux>sd*n_sd)=NaN;
%  plot(SL)
% n_std=[5 5 4 4 4 3];
% k=[30*24*60 5*24*60 3*60 60 30 5];
n_std=[5 5 4 4 3];
k=[3*60 3*60 60 30 5];
% eliminam spikes
for n=1:length(k)
    SL=spikes_v2(SL,k(n),n_std(n));
%      plot(SL)
end
%  legend('0','1','2','3','4','5','6','7')
end