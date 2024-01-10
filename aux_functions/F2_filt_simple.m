% F2_filt:
% Filtre + selecció d'episodis
function [y]=F2_filt_simple(time,y,f_altaf,f_baixaf)
%--------------------------------------------------------------------------
%--- Default parametes
%--------------- parametres filtratge -------------------------------------
% f_altaf=2; %--- minuts
% f_baixaf=3*60; %--- minuts
%--------------------------------------------------------------------------
%--- tractament de frequencia ---------------------------------------------
%--------------------------------------------------------------------------
order=4;
%--- llevam s'alta freqüència
y=butter_filter_v1(y,mode(diff(time))*24*60,f_altaf,order,'low');

%--- llevam sa baixa freqüència
y=butter_filter_v1(y,mode(diff(time))*24*60,f_baixaf,order,'high');

end