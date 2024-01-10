addpath './aux_functions/'
%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
folder=dir('.\meteotsunamis\meteo*');
files={folder(:).name}';
dates=cellfun(@(x) strrep(x,'meteotsunami_',''),files,'UniformOutput',0);


%%
factor=0.5;
radi=20;
sig=0.1;


for nd=1:length(dates);
data=dates{nd}
folder=['.\meteotsunamis\meteotsunami_',data];


load([folder,'/atm_info.mat'],'lat','lon','Patm','mareografs','time')


%%
estacions=mareografs;
dt=mode(diff(time)*24*60*60);
%%
%--- time intervals
tini=datenum(2021,6,18,16,0,0); tfi=datenum(2021,6,20,4,0,0);
tini=time(1); tfi=time(end);
% tini=datenum(2021,5,23,7,0,0); tfi=datenum(2021,5,24,3,0,0);
DT=30; 
lw=2*60;
up=[tini+lw/24/60:DT/24/60:tfi];
down=[tini:DT/24/60:tfi-lw/24/60];



%--- arreglam le sposicions
aux=isnan(lon)|isnan(lat); 
lat(aux)=[]; lon(aux)=[]; Patm(:,aux)=[]; estacions(aux)=[];

%-------------------------------------------------------------------------- 
%%
%--- definim els vetor per guardar les velocitats
vs=ones(length(up),length(estacions));
phis=vs;

%--- Mesurem lesdistàncies
Ne=length(lat);
Dist=zeros(Ne,Ne);
for ne=1:Ne
Dist(ne,:)=deg2km(distance(lat(ne),lon(ne),lat,lon));
end

%%


% limP=1./(0.8*[1:10]/60);
% limP=[120 70 35 15 2];
% limP=flip(2*2.^([0:5]*7/5));
% periods=movmedian(limP,2);
% periods=periods(2:end);

%--- Col·locació del figures
lm=0.07;
bm=0.07;
sp=0.07;
wid=(1-2*lm-0.03)/2;
hi=(1-3*sp-bm)/3;
%%
tic

df=3;
P=Patm;


corr_mat=zeros(length(estacions),length(estacions),length(up));
mlag_mat=zeros(length(estacions),length(estacions),length(up));

    for ne=1:length(estacions)
    P(:,ne)=F2_filt_simple(time,Patm(:,ne),3,120);
    end
    disp('Filtrat')
%     P=P_filt;
    %--- Eliminem les estacions que no no arriben acert valor de la
    % variança. 
    stds=std(P,1,'omitnan');
    
for ni=1:length(up)
    aux=time>=down(ni) & time<up(ni);
    fprintf('progress: %2.2f \n',ni/length(up)*100)
    P_t=P(aux,:);
    
    %--- VAriança d'aquest tram
    sd=std(P_t,1,'omitnan');
    aux=stds*factor>sd;
    
    P_t(:,aux)=NaN;
    
[est_corr_max,est_corr_mlag]=F3_correlation_v6(P_t,Dist,radi);
corr_mat(:,:,ni)=est_corr_max;
mlag_mat(:,:,ni)=est_corr_mlag;

v=zeros(length(lat),2);
m_lon=zeros(length(lat),1);
m_lat=zeros(length(lat),1);
for ne=1:length(lat)
    d=Dist(:,ne); aux=d<=radi;
    if length(find(aux))<3
        v(ne,:)=NaN;
        m_lon(ne)=NaN; m_lat(ne)=NaN;
    else
        m_lon(ne)=mean(lon(aux)); m_lat(ne)=mean(lat(aux));
        v(ne,:)=orlic_v_v1(lon(aux),lat(aux),est_corr_max(aux,aux),est_corr_mlag(aux,aux),sig,'dtime',dt,'corr_th',0.6);
    end
end

vs(ni,:)=v(:,1)';
phis(ni,:)=v(:,2)';

end
t=toc
vvv=vs;
phipp=phis;

tt=(up+down)/2;

save([folder,'/prop_vel_corr_allfreq_',data,'.mat'],'estacions','lon','lat','vvv','phipp','up','down','tt','lw','DT')
end

return
