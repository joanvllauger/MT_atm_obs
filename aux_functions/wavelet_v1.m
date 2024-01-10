% Wavelet spectrum
% Joan Villalonga 30/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [PW,period,coi,sig]=wavelet_v1(y,dt,varargin) 
%
% INPUTS:
%       y: serie temporal
%       dt: periode de mostreig
%
%   opcionals
%
%       ConfLevel: conficence level. Default 95%
%
%       lag1: fisrt autocorrelation lag. If lag1=0 it si computed here.
%       Default lag1=0.
%
%       DJ: the spacing between discrete scales. Default is 0.25.
%         A smaller # will give better scale resolution, but be slower to plot.
%
%       MotherWl: he mother wavelet function. The choices are 'MORLET', 'PAUL', or 'DOG'
%
%       PARAM = the mother wavelet parameter.
%             For 'MORLET' this is k0 (wavenumber), default is 6.
%             For 'PAUL' this is m (order), default is 4.
%             For 'DOG' this is m (m-th derivative), default is 2.
%
%       PAD = if set to 1 (default is 1), pad time series with enough zeroes to get
%          N up to the next higher power of 2. This prevents wraparound
%          from the end of the time series to the beginning, and also
%          speeds up the FFT's used to do the wavelet transform.
%          This will not eliminate all edge effects (see COI below).
%
%       S0 = the smallest scale of the wavelet.  Default is 2*DT.
% 
%       J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%         to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
% 
%
% OUTPUTS:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PW,period,coi,sig]=wavelet_v1(y,dt,varargin)
%--- Dafault parameters
MotherWL='Morlet'; Dj=1/16; PARAM=12; ConfLevel=0.95; PAD=1; lag1=0;
S0=2*dt; 
Tmax=360; % minuts
N=length(y);
J1 =log2(N*dt/S0)/Dj;
NA=1;
S0=-1; J1=-1;

%--- Assiganam les variables opcionals 
if (mod(length(varargin),2)==1), error('Number of imputs has to be even'); 
else
%--- agafam els imparells del varargin que són els que duran els noms
in=1:2:length(varargin); noms=varargin(in);
%--- agafam els parells que son els que tenen els valors dels par?metres
in=2:2:length(varargin); valors=varargin(in); clear in;
end

%--- Comprovam que tots els par?metres siguin correctes.
PosParam={'MotherWL','Dj','PARAM','ConfLevel','J1','S0','lag1','NA','Tmax'};
%--- Comprovam que cada valor dels par?metres s?n del tipo que toca
PosParamType={'char','double','double','double','double','double','double','double','double'};
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
if J1==-1
    s0=2*dt; 
J1=log2((Tmax*(PARAM+sqrt(2+PARAM^2)))/(4*pi*s0))/Dj;
J1=ceil(J1);
end

%--------------------------------------------------------------------------
%- Wavelet transform de les wavelets
[wt,period,scale,coi] = wavelet_trans(y,dt,PAD,Dj,S0,J1,MotherWL,PARAM);
PW=(abs(wt).^2);



if lag1==0
[lag, corr]=correlacion_lag(y,y,1);
lag1=corr(end);
end
if NA==1 % alerta si es canvia NA perquè hi pot haver problemes
dof=-1;
else
dof=NA;    
end

%- Wavelet significance 
[signif,fft_theor] = wave_signif(y,dt,scale,0,lag1,ConfLevel,dof,MotherWL,PARAM);
sig = (signif')*(ones(1,N));  % expand signif --> (J+1)x(N) array
sig = PW./sig;  
if NA==1
else
PW=movmean(PW,NA,'omitnan');
end
end