% Funció del càlcul de velocitat
% Joan Villalonga
% 02/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [v]=orlic_v_v1(lon,lat,est_corr_max,est_corr_mlag,sig_c,varargin) 
%
% INPUTS:
%   lon: vector amb les longituds de totes les estacions
%   lat: vector amb les latituds de totes les estacions
%   est_corr_max: matriu amb la correlació màxima entre estacions (size=N_e x N_e)
%   est_corr_mlag: matriu amb el lag on es troba el màxim de correlació entre estacions (size=N_e x N_e)
%   sig_c: llindar de correlació per el que el pes és l'1%
%
% opcional:
%   dtime: diferencial de temps en segons per escalar el lag a segons 
%
% OUTPUTS:
%   v: 
%       v(1)--> c: modul de la velocitat
%       v(2)--> gamma: angle del vector velocitat respecte als paral·lels
%--------------------------------------------------------------------------
function [c]=orlic_v_v1(lon,lat,est_corr_max,est_corr_mlag,sig,varargin)
%--- Definim variables globals 
global corr_m gamma_i Dt di sig_c
sig_c=sig;
corr_th=0.6;
%--- Parametres editables
dtime=60; %segons

%--- Assiganam les variables opcionals 
if (mod(length(varargin),2)==1), error('Number of imputs has to be even'); 
else
%--- agafam els imparells del varargin que sÛn els que duran els noms
in=1:2:length(varargin); noms=varargin(in);
%--- agafam els parells que son els que tenen els valors dels par?metres
in=2:2:length(varargin); valors=varargin(in); clear in;
end

%--- Comprovam que tots els par?metres siguin correctes.
PosParam={'dtime','corr_th'};
%--- Comprovam que cada valor dels par?metres s?n del tipo que toca
PosParamType={'double','double'};
for i=1:length(noms)
    if contains(noms{i},PosParam)==0, error(['The parameter ',noms{i},' is not valid.']); 
    else
        k=find(strcmp(noms{i},PosParam));
        if strcmp(class(valors{i}),PosParamType{k})==0 
        error(['The value of ',noms{i},' must be ',PosParamType{k},' class, not ',class(valors{i}),' class']); 
        else
        eval([noms{i},'=valors{i};']);
        end
    end
end
lon=abs(lon);

%--- Eliminam estacions que no sabem on soon
aux=not(isnan(lon) | isnan(lat)); 
lon=lon(aux); lat=lat(aux); est_corr_max=est_corr_max(aux,aux); est_corr_mlag=est_corr_mlag(aux,aux); 
%--- Fem la matriu de direfències de longitud i latitud
dlat=lat'-lat;
dlon=lon'-lon;
%--- Fem una matriu amb el cosinus de totes les latituds
cos_phi=cos(pi/180*ones(length(lat),length(lat)).*lat');
%--- Calculem la matriu de distàncies entre estacions
Dist=zeros(length(lat),length(lat));
for ne=1:length(lat)
Dist(ne,:)=deg2km(distance(lat(ne),lon(ne),lat,lon));
end

%--- Utilitzem tant sols els punts on el lag és positiu
aux=find(est_corr_mlag>0);
di=Dist(aux); corr_m=est_corr_max(aux); corr_lag=est_corr_mlag(aux);
dlon=dlon(aux); dlat=dlat(aux); cos_phi=cos_phi(aux);

%--- Si tenim nans els eleminam
aux=not(isnan(corr_m)|isnan(corr_lag)|corr_m<corr_th);
di=di(aux); corr_m=corr_m(aux); corr_lag=corr_lag(aux);
dlon=dlon(aux); dlat=dlat(aux); cos_phi=cos_phi(aux);
% 
% di
% corr_lag
% corr_m


if length(corr_m)<3
    c=[NaN NaN];
else

%--- Analisi per quadrants i obtenció de gamma
gamma_i=zeros(length(dlon),1);
for k=1:length(dlon)
    dlat_i=dlat(k);
    dlon_i=dlon(k);
%- 1rQ
if dlat_i>0 & dlon_i>0
    gamma_i(k)=atan(dlat_i/(dlon_i*cos_phi(k)));
%- 2nQ    
elseif dlat_i>0 & dlon_i<0
    gamma_i(k)=pi+atan(dlat_i/(dlon_i*cos_phi(k)));
%- 3rQ    
elseif dlat_i<0 & dlon_i<0
    gamma_i(k)=pi+atan(dlat_i/(dlon_i*cos_phi(k)));
%- 4tQ     
elseif dlat_i<0 & dlon_i>0
    gamma_i(k)=2*pi+atan(dlat_i/(dlon_i*cos_phi(k)));
end
end
    
%--- Passem a les unitats adequades
Dt=corr_lag*dtime; % passam a segons
di=di*1000; % passam a metres

x0=[10;pi/4];
ops=optimset('MaxFunEvals',1000000,'MaxIter',1000000);
c=fminsearch(@dt_error_wcorr,x0,ops);
c(2)=mod(c(2),2*pi);
end
end