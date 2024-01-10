function [t,wh,t_wh,wh_raw]=wave_height_v1(time,y,varargin)
%--------------------------------------------------------------------------
%--- Default parametes
%--------------- parametres filtratge -------------------------------------
f_filt=1; %--- Flag flitre 0->no filtre 1-> filtre
f_altaf=2; %--- minuts
f_baixaf=2*60; %--- minuts

%--------------- parametres envolvent -------------------------------------
finestra=60; %---minuts

%--- Assiganam les variables opcionals 
if (mod(length(varargin),2)==1), error('Number of imputs has to be even'); 
else
%--- agafam els imparells del varargin que són els que duran els noms
in=1:2:length(varargin); noms=varargin(in);
%--- agafam els parells que son els que tenen els valors dels par?metres
in=2:2:length(varargin); valors=varargin(in); clear in;
end

%--- Comprovam que tots els par?metres siguin correctes.
PosParam={'f_filt','f_altaf','f_baixaf','finestra'};
%--- Comprovam que cada valor dels par?metres s?n del tipo que toca
PosParamType={'double','double','double','double'};
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

x=y;
%--- llevam nan inicials y finals
in=find(isnan(y)==0,1,'first'); fi=find(isnan(y)==0,1,'last');
time=time(in:fi); y=y(in:fi);
nans=isnan(y);
%--- cosa extra
% limit=10;
% in=find(cumsum(not(isnan(y)))>limit,1,'first');
% time=time(in:end); y=y(in:end);

%--- interpolam linealment els forats
y=interp1(time(not(isnan(y))),y(not(isnan(y))),time);

%--------------------------------------------------------------------------
%--- tractament de frequencia ---------------------------------------------
%--------------------------------------------------------------------------
if logical(f_filt)
%--- llevam s'alta freqÃ¼encia
y=butter_filter_v1(y,mode(diff(time))*24*60,f_altaf,2,'low');

%--- llevam sa baixa freqÃ¼encia
y=butter_filter_v1(y,mode(diff(time))*24*60,f_baixaf,2,'high');
end
y(nans)=NaN;
% figure
% plot(datetime(time,'ConvertFrom','datenum'),y)

[pks,locs]=findpeaks(abs(y));
% figure
% plot(time,abs(y))
% hold on
% plot(time,y)
% plot(time(locs),y(locs),'o')


wh_raw=diff(y(locs));
t_wh=time(locs(2:end));

%%
%--- definim una finestra per cercar mÃ xims
dw=finestra;
ww=dw/(2*24*60);
t=[time(1)+ww:2*ww:time(end)]';
t=round(t*24)/24;
wh=zeros(length(t),1);
for nt=1:length(t)
    aux=(t_wh>=t(nt)-ww & t_wh<t(nt)+ww);
        whm=max(abs(wh_raw(aux)));
        if isempty(whm)
        wh(nt)=NaN;
        else
             wh(nt)=whm;
        end

end

end
