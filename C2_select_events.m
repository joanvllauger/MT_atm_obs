addpath './aux_functions/'
% codi per trobar l'inici i final de tots els esdeveniments de risssaga que
% tenc
clear all

% 1--- Legim les dades de nivell del mar de Ciutadalla
load('./data/ciutadella_SL_AtmPres.mat','time','SL')

% 2-- We cut the time series for the period of study
SL=SL/100;

[t,wh,t_wh,wh_raw]=wave_height_v1(time,SL);

%%
durada_events=1; % dia

%%
th_max=0.6;
aux=wh>th_max;
wh_events=wh(aux);
t_events=t(aux);

th_max=0.6;
aux=wh>th_max;
wh_events=wh(aux);
t_events=t(aux);
c=1;
event=[];
wh_event=[];
while c<=length(t_events)
    aux=find(abs(t_events(c)-t_events(c:end))<durada_events)+c-1;
    [mm,id]=max(wh_events(aux));
    wh_event=[wh_event;mm];
    event=[event; t_events(aux(id))];
    c=c+length(aux);
end

%%
th_min=prctile(wh,90);
t_in=zeros(length(wh_event),1);
t_out=t_in;
for n=1:length(wh_event)
    aux=find(t==event(n));
    in=find(wh(1:aux)<th_min,1,'last');
    taux=t(1:aux);
    t_in(n)=taux(in);
    out=find(wh(aux:end)<th_min,1,'first');
    taux=t(aux:end);
    t_out(n)=taux(out);
    
end

%%
c=1;
event_new=[];
wh_event_new=[];
tfis=[];
tins=[];
while c<=length(event)
    aux=find(round((t_in(c:end)-t_out(c))*24)<0)+c-1;
    [mm,id]=max(wh_event(aux));
    wh_event_new=[wh_event_new;mm];
    event_new=[event_new; event(aux(id))];
    tfis=[tfis; max(t_out(aux))];
    tins=[tins; min(t_in(aux))];
    
    if isempty(aux)
        c=c+1;
    else
    c=c+length(aux);
    end
end

aux=event_new>=datenum(2021,1,1) & event_new<datenum(2022,1,1);
T=table(datetime(event_new(aux),'ConvertFrom','datenum'),datetime(tins(aux),'ConvertFrom','datenum'),datetime(tfis(aux),'ConvertFrom','datenum'),wh_event_new(aux));
writetable(T,'List_of_events.xlsx')

T=table(datetime(event_new,'ConvertFrom','datenum'),datetime(tins,'ConvertFrom','datenum'),datetime(tfis,'ConvertFrom','datenum'),wh_event_new);
writetable(T,'List_of_events_long.xlsx')


