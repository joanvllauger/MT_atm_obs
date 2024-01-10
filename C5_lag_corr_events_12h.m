close all; clear all;
addpath './aux_functions/'
% dibuixam les correlacions 
load(['./data/Atm_pres_all.mat'])

aux=time>=datenum(2021,4,1) & time<datenum(2021,11,1);
time=time(aux);
Patm=Patm(aux,:);


%--- Importam les rissagues que tenim controlades
T=readtable('List_of_events.xlsx');
times=datenum(table2array(T(:,1)));

%--- Carregem les longitud i latitud bones
% T=readtable('Patm_stations_list.xlsx');
% lon=table2array(T(:,2));
% lat=table2array(T(:,3));

wint=12;

% 1- Eliminar spikes

% 2- Eliminam les mesures que son 0
aux=Patm==0; Patm(aux)=NaN;

%--- Eliminam els salts massa grans
for ne=1:length(mareografs)
    dd=diff(Patm(:,ne));
    aux=find(abs(dd)>1);
    for nf=1:length(aux)
        Patm(aux(nf),ne)=NaN;
        Patm(aux(nf)+1,ne)=NaN;
    end
end
% 3- Eliminar ourliers
finestra=10; % min
nanlen=30; % min
%--- Calculam la desviació std per la finestra
stds=movstd(Patm,finestra,0,1,'omitnan');
%--- definim un threshold en std
th_std=0.0000004;
%--- Eliminem els punts amb variança molt petita
aux=stds<th_std; Patm(aux)=NaN;

%--- Eliminam punts aillats
disp('eliminem dades aillades')
for ne=1:length(mareografs)
    nonan=find(not(isnan(Patm(:,ne))));
    nan=find((isnan(Patm(:,ne))));
    f_nonan=nonan(find(not(diff(nonan)==1))+1);
    f_nan=nan(find(not(diff(nan)==1))+1);
    if isnan(Patm(1,ne))
        f_nonan=[nonan(1);f_nonan];
    else
        f_nan=[nan(1);f_nan];
        f_nonan=[1;f_nonan];
    end
    nn=min(length(f_nan),length(f_nonan));
    len=f_nan(1:nn)-f_nonan(1:nn);    
    for nl=1:length(len)
        if len(nl)<nanlen
            Patm(f_nonan(nl):f_nan(nl),ne)=NaN;
        end
    end

end

%% Filtram la serie

%--- Definiem la variable
mat=Patm; %--- temps x espacions
aux=isnan(mat);
%--- Eliminem la mitjana 
for ne=1:length(mareografs)
    mat(:,ne)=mat(:,ne)-mean(mat(:,ne),'omitnan');
    if isnan(mat(1,ne))
        mat(1,ne)=0;
    end
    if isnan(mat(end,ne))
        mat(end,ne)=0;
    end
end
%--- Omplim forats
mida_gap_max=12*60*1e10; % min
for ne=1:length(mareografs)
mat(:,ne)=rellenar_huecos(time,mat(:,ne),mida_gap_max/24/60,'lineal');
end
%--- Filtram
mat_0=mat;
for ne=1:length(mareografs)
 %    mat_0(:,ne)=F2_filt_simple(time,mat_0(:,ne),2,2*60);
    mat_0(:,ne)=butter_filter_v1(mat_0(:,ne),1,2*60,2,'high');
end
mat=mat_0; mat(aux)=NaN;
P_filt=mat;

disp('Proces d''eliminació de nans acabat')
%% Eliminem els lats en els forats
% P_filt_old=P_filt;
th_diff=1; % hPa;
lennan=1.5*60; % minutes
%--- Eliminam punts aillats
disp('Eliminem les dades dels gaps on hi ha salts')
for ne=1:length(mareografs)
    nonan=find(not(isnan(P_filt(:,ne))));
    nan=find((isnan(P_filt(:,ne))));
    f_nonan=nonan(find(not(diff(nonan)==1))+1);
    f_nan=nan(find(not(diff(nan)==1))+1);
    if not(isnan(P_filt(1,ne)))
        f_nan=[nan(1);f_nan];
    end
    nn=min(length(f_nan),length(f_nonan));
    dd=abs(P_filt(f_nan(1:nn)-1,ne)-P_filt(f_nonan(1:nn),ne));    
    for nl=1:length(dd)
        if dd(nl)>th_diff
            if f_nan(nl)-lennan<1
                P_filt(1:f_nonan(nl)+lennan,ne)=NaN;
            elseif f_nonan(nl)+lennan>length(P_filt(:,ne))
                P_filt(f_nan(nl)-lennan:end,ne)=NaN;
            else
                P_filt(f_nan(nl)-lennan:f_nonan(nl)+lennan,ne)=NaN;
            end
        end
    end
end
disp('Serie Filtrada areglada')


%% Calcul de les matrius de correlació

%---- dividim la matriu de pressió en els trossos necessaris
P_t=cell(length(times),1);
P_nofilt=cell(length(times),1);
for nt=1:length(times)
    aux=time>=times(nt)-wint/2/24&time<times(nt)+wint/2/24;
    P_t{nt}=P_filt(aux,:);
    P_nofilt{nt}=Patm(aux,:);
    
end
corr_max_cell=cell(length(times),1); 
corr_mlag_cell=cell(length(times),1); 
pvalue_cell=cell(length(times),1); 

%--- Feim la matriu de distàncies
Ne=length(mareografs);
Dist=zeros(Ne,Ne);
for ne=1:Ne
Dist(ne,:)=deg2km(distance(lat(ne),lon(ne),lat,lon));
end
%%
disp('entram al parfor')
parfor nt=1:length(times)

tic
P=P_t{nt};

[corr_max_mat,corr_mlag_mat,pvalue_mat]=F3_correlation_v4(P,Dist);

corr_max_cell{nt}=corr_max_mat;
corr_mlag_cell{nt}=corr_mlag_mat; 
pvalue_cell{nt}=pvalue_mat;
toc

disp(['Iteration: ',num2str(nt/length(times))]) 

end
save(['./data/corr_rissagues_1min_allfreq_12h.mat'],'corr_max_cell','pvalue_cell','corr_mlag_cell','times','lon','lat','mareografs','time','P_filt','P_t','P_nofilt','Patm','-v7.3')

return
%%
load(['./data/corr_rissagues_1min_allfreq_12h.mat'],'corr_max_cell','pvalue_cell','corr_mlag_cell','times','lon','lat','mareografs','time','P_filt','P_t','P_nofilt','Patm')

%%
n=3;
for n=1:length(mareografs)
figure
plot(P_t{1}(:,n))
end




