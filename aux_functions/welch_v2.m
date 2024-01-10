% Welch-Berlett Perididogram
% Joan Villalonga 13/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2: ens dona la significan�a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUTS:
%       y: series temporal
%       WindowLength: longitud de la finestra apicada
%       dt: periode de mostreig
%
%   opcionals
%       PercentOverlap: % de superposici� de les finestres. Entre 0 i 100.
%       Confidence: tant per 1 de confian�a. 
%       WindowType: tipus de finestra, es recomana 'hamming'o 'kaiser' (inlcoure el nomque te la funci� de Matlab)
%       Beta: (nom�s si s'indica la finsetra de keiser) amplada de la finestra de Kaiser.
%   valors per defecte -->
%          PercentOverlap=50; Confidence=0.95; WindowType='kaiser'; Beta=10;
% OUTPUTS:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,f,conf,DOF]=welch_v2(y,WindowLength,dt,varargin)

%--- Variables obcional per defecte
PercentOverlap=50; Confidence=0.95; WindowType='kaiser'; Beta=10;
WindowType='hamming';
%--- Assiganam les variables opcionals 
if (mod(length(varargin),2)==1), error('Number of imputs has to be even'); 
else
%--- agafam els imparells del varargin que s�n els que duran els noms
in=1:2:length(varargin); noms=varargin(in);
%--- agafam els parells que son els que tenen els valors dels par?metres
in=2:2:length(varargin); valors=varargin(in); clear in;
end

%--- Comprovam que tots els par?metres siguin correctes.
PosParam={'PercentOverlap','Confidence','WindowType','Beta'};
%--- Comprovam que cada valor dels par?metres s?n del tipo que toca
PosParamType={'double','double','char','double'};
for i=1:length(noms)
    if contains(noms{i},PosParam)==0, error(['The parameter ',noms{i},' is not valid.']); 
    elseif strcmp(class(valors{i}),PosParamType{i})==0 
        error(['The value of ',noms{i},' must be ',PosParamType{i},' class, not ',class(valors{i}),' class']); 
    else
        eval([noms{i},'=valors{i};']);
    end
end

%--------------------------------------------------------------------------

%--- Feim la finestra
if strcmp(WindowType,'kaiser')
    window=kaiser(WindowLength,Beta);
else
    eval(['window=',WindowType,'(WindowLength);']);
end

%--------------------------------------------------------------------------
%--- Optenim freq��ncia de mostreig
Fs=1/dt;

%--- Optenim noverlap
noverlap=round(PercentOverlap/100*WindowLength);

%--- Calcul del perdidiograma
[P,f,conf] = pwelch(y,window,noverlap,[],Fs);%,'ConfidenceLevel',Confidence);

%--------------------------------------------------------------------------
%--- Calculem nfft
M=2*(length(f)-1);

%--- Calcul de la Significan�a
N=length(y);
%-- Graus de llibertat
DOF=2*(floor(N/M)+floor((N-M/2)/M));

end

