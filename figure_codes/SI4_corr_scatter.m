%--- Ubicació de les dades a obrir
addpath '../aux_functions/'

path_mat='..\meteotsunamis\meteotsunami_';

%--- Inicialitzem un interval de temps que ens interesa
tini=datenum(2021,5,23);
load(['..\meteotsunamis\meteotsunami_',datestr(tini,'yyyymmdd'),'/prop_vel_corr_',datestr(tini,'yyyymmdd'),'.mat'],'periods','limP');
load(['..\meteotsunamis\meteotsunami_',datestr(tini,'yyyymmdd'),'/atm_info.mat'],'time','Patm','lon','lat','mareografs');
estacions=mareografs;

DDist=zeros(length(lon),length(lon));
for ne=1:length(lat)
DDist(ne,:)=deg2km(distance(lat(ne),lon(ne),lat,lon));
end


k=1;

FS=10;
MS=8;

nh=2;
nv=2;
lm=0.07;
rm=0.02;
bm=0.08;
tm=0.045;
hs=0.03;
vs=tm;

hi=(1-(nv-1)*vs-tm-bm)/nv;
wi=(1-(nh-1)*hs-lm-rm)/nh;
lletres=flip({'a','b','c','d'});
figure('units','centimeters','Position',[2 2 17 13])
for np=1:length(periods)-1
    nf=np;
    freq=periods(np);
    for ne=1:length(estacions)
    P(:,ne)=F2_filt_simple(time,Patm(:,ne),limP(np+1),limP(np));
    end
    disp('Filtrat')

    
    
[est_corr_max,est_corr_mlag]=F3_correlation_v6(P,DDist,200);

h=double(mod(k,2));
v=1;
if k>2
    v=2;
end


% aux=not(isnan(lon) | isnan(lat)); 
% Dist=Dist(aux,aux);
% lon=lon(aux); lat=lat(aux); est_corr_max=est_corr_max(aux,aux); est_corr_mlag=est_corr_mlag(aux,aux); 


%--------------------------------------------------------------------------
%--- Angles
%--- Eliminam estacions que no sabem on soon
aux=not(isnan(lon) | isnan(lat)); 
Dist=DDist(aux,aux);
lon=lon(aux); lat=lat(aux); est_corr_max=est_corr_max(aux,aux); est_corr_mlag=est_corr_mlag(aux,aux); 
%--- Fem la matriu de direfències de longitud i latitud
dlat=lat'-lat;
dlon=lon'-lon;
%--- Fem una matriu amb el cosinus de totes les latituds
cos_phi=cos(pi/180*ones(length(lat),length(lat)).*lat');

%--- Analisi per quadrants i obtenció de gamma
gamma_i=zeros(size(dlon));
for nn=1:length(dlon)
    for kk=1:length(dlon)
    dlat_i=dlat(nn,kk);
    dlon_i=dlon(nn,kk);
        %- 1rQ
        if dlat_i>0 & dlon_i>0
            gamma_i(nn,kk)=atan(dlat_i/(dlon_i*cos_phi(nn,kk)));
        %- 2nQ    
        elseif dlat_i>0 & dlon_i<0
            gamma_i(nn,kk)=pi+atan(dlat_i/(dlon_i*cos_phi(nn,kk)));
        %- 3rQ    
        elseif dlat_i<0 & dlon_i<0
            gamma_i(nn,kk)=pi+atan(dlat_i/(dlon_i*cos_phi(nn,kk)));
        %- 4tQ     
        elseif dlat_i<0 & dlon_i>0
            gamma_i(nn,kk)=2*pi+atan(dlat_i/(dlon_i*cos_phi(nn,kk)));
        end
    end
end

%--------------------------------------------------------------------------

glim1=(270-210)*pi/180;
glim2=(270-260)*pi/180;

glim3=(450-210)*pi/180;
glim4=(450-260)*pi/180;

auxg=gamma_i>=glim2 & gamma_i<=glim1 | gamma_i>=glim4 & gamma_i<=glim3;


[N,EDGES,BIN] = histcounts(Dist,30);

axes('Position',[lm+(hs+wi)*h bm+(v-1)*(hi+vs) wi hi])
r=plot(Dist(not(auxg)),est_corr_max(not(auxg)),'.r','MarkerSize',MS);
hold on
dd=movsum(EDGES,2)/2;
dd=dd(1:end);
my=ones(length(dd),1);
p25=ones(length(dd),1);
p75=ones(length(dd),1);
for nb=2:length(dd)
    y=est_corr_max(BIN==nb-1);
    my(nb)=median(y);
    p25(nb)=prctile(y,25);
    p75(nb)=prctile(y,75);
    
end
b=plot(Dist(auxg),est_corr_max(auxg),'.b','MarkerSize',MS);

    
   
% b=plot(dd,my,'b','Linewidth',3);
% plot(dd,p25,'--r','Linewidth',3);
% plot(dd,p75,'-.r','Linewidth',3);
ax=gca;
ax.FontSize=FS;
xlim([0 200]);


if h
ylabel('')
ax.YTickLabel='';
else
    ylabel('Correlation coefficient') 
end
ax=gca;
ticks=[0:25:200];
ax.XTick=ticks;
box=[0.06 0.03];
grid on

ylim([0 1])
ax.YTick=[0:0.1:1];

if k<3
    xlabel('Distance [km]')
else
    ax.XTickLabel='';
end

ax.YMinorTick='off';
ax.YMinorGrid='off';
ax.XMinorTick='off';
ax.XMinorGrid='off';
ax.GridAlpha = 0.8 ;
ax.MinorGridColor = 'k';
ax.MinorGridAlpha = 0.8 ;

if k==4
    legend([b,r],'Direction of propagation','Other directions')
end

title(['(',lletres{nf},')  ', num2str(limP(nf+1),'%3.0f'),' to ',num2str(limP(nf),'%3.0f'),' minutes period'])
k=k+1;
end

% saveas(gcf,['./TFM_figs/Patm_dist_corr_',num2str(freq,'%3.0f'),'_',datestr(tini,'mmdd')],'png')
print(gcf,'-djpeg','-r500',['SI4'])
