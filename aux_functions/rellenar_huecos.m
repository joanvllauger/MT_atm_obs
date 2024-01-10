function [ynew]=rellenar_huecos(time,y,tmax_gap,tipo)
%--- encuentro los huecos y su longitud
aux=isnan(y);
if sum(aux)==0
    ynew=y;
else
ingap=find(diff(aux)==1)+1;
figap=find(diff(aux)==-1);
if not(isempty(ingap)) & isempty(figap)
    figap=length(y);
end
if isempty(ingap) & not(isempty(figap))
    ingap=1;
end
if (ingap(1)>figap(1)), ingap=cat(1,1,ingap); end
if (ingap(end)>figap(end)), figap=cat(1,figap,length(y)); end

%--- numero inicial de huecos
NHin=length(ingap);

%--- haremos uno de los dos tipos de interpolacion
switch tipo
    case 'lineal'
        ynew=interp1(time(not(isnan(y))),y(not(isnan(y))),time,'linear',0);
    case 'autorregresivo'
        npuntos=round(100*tmax_gap/mode(diff(time)));
        ynew=fillgaps(y,npuntos);
    case 'pchip'
        ynew=interp1(time(not(isnan(y))),y(not(isnan(y))),time,'pchip',0);
end

%--- ahora tengo que poner nan a los huecos más largos de lo que quiero
%rellenar
NHfi=NHin;
NP=round(tmax_gap/mode(diff(time)));
for nh=1:length(ingap)
    if ((figap(nh)-ingap(nh)+1)>NP) %--- el hueco se queda
        ynew(ingap(nh):figap(nh))=NaN;
    
    else %--- el hueco se queda, tendré un hueco menos
        NHfi=NHfi-1;
    end
end

% disp(['# inicial gaps = ',num2str(NHin),'; # final gaps = ',num2str(NHfi),...
%     '; % gaps eliminado = ',num2str(100*(1-NHfi/NHin),'%2.2f'),'%']);

end
end

