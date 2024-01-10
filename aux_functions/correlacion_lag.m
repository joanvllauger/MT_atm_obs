% [lag, corr]=correlacion_lag(x,y,maxlag)
% Funció de correlació entre dues variables feta per Àngle Amores
function [lag, corr]=correlacion_lag(x,y,maxlag)
%% Lo transformamos todo a columnas
if (size(x,2)~=1), x=x'; end; if (size(y,2)~=1), y=y'; end
%% Le quitamos las medias
x=x-nanmean(x); y=y-nanmean(y);
%% nos creamos la variables iniciales
lag=(-maxlag:1:maxlag)'; corr=zeros(size(lag)); sig=zeros(size(lag));

%% Correlación una parte
for i=1:(length(lag)-1)/2
    auxX=x(abs(lag(i))+1:end); auxY=y(1:end-abs(lag(i)));
    %correlación
    aux=nancov([auxX,auxY]); %corr(i)=aux(2,1)./sqrt(sigmax*sigmay); clear aux;
    corr(i)=aux(2,1)./sqrt(aux(1,1)*aux(2,2)); clear aux;
end
%la otra
for i=(length(lag)+1)/2:length(lag)
    auxX=x(1:end-abs(lag(i))); auxY=y(abs(lag(i))+1:end);
    %correlación
    aux=nancov([auxX,auxY]); %corr(i)=aux(2,1)./sqrt(sigmax*sigmay); clear aux;
    corr(i)=aux(2,1)./sqrt(aux(1,1)*aux(2,2)); clear aux;
end
