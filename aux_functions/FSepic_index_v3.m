% Funció amb el codi den Biel per Calcular l'Index Sepic 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   time_era: vector de temps de les dades del reanalisi
%   var_era:  matriu [temps x vars]
%   t_wh:     vector de temps de les d'alçades d'ona
%   wh:       vector d'alçades d'ona
%
% OUTPUT
%   wh_filt: alçades d'ona filtrades
%   S_ind: Sepic index
%   parameters: parametres del fitting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wh_filt,S_ind,SMOS_ind,parameters]=FSepic_index_v3(time_era,var_era,t_wh,wh,varargin)
%--- parameters
windays=5; % in days
smoothmax=0; % If instead of using mean in 5 days we use 95 quantile value - NOT WORTH IT
smoothmax_era=0;
%--- Assiganam les variables opcionals 
if (mod(length(varargin),2)==1), error('Number of imputs has to be even'); 
else
%--- agafam els imparells del varargin que són els que duran els noms
in=1:2:length(varargin); noms=varargin(in);
%--- agafam els parells que son els que tenen els valors dels par?metres
in=2:2:length(varargin); valors=varargin(in); clear in;
end

%--- Comprovam que tots els par?metres siguin correctes.
PosParam={'windays','smoothmax'};
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


% 3day filtering
dt=median(diff(time_era));
if windays==1/24
    SMOVAR=var_era;
    VRED=NaN+zeros(length(t_wh),size(var_era,2));
    for nk=1:size(var_era,2)
    VRED(:,nk)=interp1(time_era,SMOVAR(:,nk),t_wh);
    end
    iok=find(~isnan(sum([VRED wh(:)],2)));
    X=VRED(iok,:);
    B=wh(iok);
    parameters=X\B;
    S_ind=VRED*parameters;
    wh_filt=wh;
    SMOS_ind=S_ind;
    disp('hola')
else
    
nwin=round(windays/dt/2)

finestra=kaiser(nwin*2+1,2);
SMOVAR=var_era.*NaN;
for nt=nwin+1:length(time_era)-nwin;
    %     if mod(nt,1000)==0;disp(sprintf('Smoothing %i of %i',nt,length(time_era)));end
%            if smoothmax==1
%             SMOVAR(nt,:)=quantile(var_era(nt-nwin:nt+nwin,:),0.95);
%            else
            SMOVAR(nt,:)=nanmean(var_era(nt-nwin:nt+nwin,:).*finestra,1);
%            end
end

% Filter the wave data
dt=median(diff(t_wh));
nwin=round(windays/dt/2);
finestra=kaiser(nwin*2+1,2);

wh_filt=t_wh.*NaN;
for nt=nwin+1:length(t_wh)-nwin;
    if smoothmax==1
%     wh_filt(nt)=quantile(wh(nt-nwin:nt+nwin),0.95);
        wh_filt(nt)=max(wh(nt-nwin:nt+nwin));
    else
        wh_filt(nt)=nanmean(wh(nt-nwin:nt+nwin).*finestra);
    end
end

% Fit the variables to the SL data
% .. new data interpolated to SL times
VRED=NaN+zeros(length(t_wh),size(var_era,2));
dw=mode(round(diff(t_wh)*24*60));
ww=dw/(2*24*60);
t=t_wh;
for nk=1:size(var_era,2)
for nt=1:length(t)
    aux=(time_era>=t(nt)-ww & time_era<t(nt)+ww);
    if smoothmax_era==1
        whm=max(abs(SMOVAR(aux,nk)));
    else
        whm=nanmean(abs(SMOVAR(aux,nk)));
    end
        if isempty(whm)
            VRED(nt,nk)=NaN;
        else
            VRED(nt,nk)=whm;
        end

end
end

% Fitting
iok=find(~isnan(sum([VRED wh_filt(:)],2)));
X=VRED(iok,:);
B=wh_filt(iok);
parameters=X\B;
S_ind=VRED*parameters;


% 3day filtering
dt=mode(diff(t_wh));
finestra=kaiser(nwin*2+1,2);

SMOS_ind=S_ind.*NaN;
for nt=nwin+1:length(t_wh)-nwin;
    %     if mod(nt,1000)==0;disp(sprintf('Smoothing %i of %i',nt,length(time_era)));end
%            if smoothmax==1
%  SMOS_ind(nt,:)=quantile(S_ind(nt-nwin:nt+nwin,:),0.95);
%            else
         SMOS_ind(nt,:)=nanmean(S_ind(nt-nwin:nt+nwin,:).*finestra,1);
%            end
end
% nans=isnan(SMOS_ind);
% y=rellenar_huecos(t_wh,SMOS_ind,t_wh(end),'lineal');
% SMOS_ind=butter_filter_v1(y,dt,15,2,'high');
% SMOS_ind(nans)=NaN;
end
end

