addpath './aux_functions/'
%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
folder=dir('.\meteotsunamis\meteo*');
files={folder(:).name}';
dates=cellfun(@(x) strrep(x,'meteotsunami_',''),files,'UniformOutput',0);
% dates(1)=[];

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
% for ne=[1 37]%:length(estacions)
% figure
% plot(time,Patm(:,ne))
% title(estacions{ne})
% end


dt=mode(diff(time)*24*60*60)
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
limP=[120 70 35 15 2];
limP=flip(2*2.^([0:5]*7/5));
periods=movmedian(limP,2);
periods=periods(2:end);

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

vvv=cell(length(periods),1);
phipp=cell(length(periods),1);

stds=zeros(length(periods),length(estacions));

corr_cell=cell(length(periods),1);
mlag_cell=cell(length(periods),1);

corr_mat=zeros(length(estacions),length(estacions),length(up));
mlag_mat=zeros(length(estacions),length(estacions),length(up));

for np=1:length(periods)
    freq=periods(np)
    for ne=1:length(estacions)
    P(:,ne)=F2_filt_simple(time,Patm(:,ne),limP(np+1),limP(np));
    end
    disp('Filtrat')
%     P=P_filt;
    %--- Eliminem les estacions que no no arriben acert valor de la
    % variança. 
    stds(np,:)=std(P,1,'omitnan');
    
    
    %{
    nn=30;
    figure
    plot(time,P(:,nn))
%     hold on
%     plot(time,P_filt(:,nn))
    title(estacions{nn})
    datetick
    P_prova=Patm(:,nn);
    %}
    
for ni=1:length(up)
    aux=time>=down(ni) & time<up(ni);
    fprintf('progress: %2.2f \n',ni/length(up)*100)
    P_t=P(aux,:);
    
    %--- VAriança d'aquest tram
    sd=std(P_t,1,'omitnan');
    aux=stds(np,:)*factor>sd;
    
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
%        v(ne,:)=triangle_v1(lon(aux),lat(aux),est_corr_max(aux,aux),est_corr_mlag(aux,aux),sig,'dtime',60);
    end
end

% for ne=1:length(lat)
%     d=Dist(:,ne); 
%     rho=est_corr_max(:,ne);
%     [dd,aux]=sort(rho,'descend');
%     
%     aux=aux(find(not(isnan(dd)),3));
%     cc=est_corr_max(aux,aux);
%     lags=est_corr_mlag(aux,aux);
% 
%     if  length(aux)<3 
%         v(ne,:)=NaN;   
%     elseif cc(1,2)+cc(2,3)-cc(1,3)>=2 || lags(1,2)==0 || lags(1,3)==0 || lags(2,3)==0 || min(cc,[],'all')<0.6 
%         v(ne,:)=NaN; 
%     else
%         [v(ne,1),lambda(ne),v(ne,2)]=triangle_v1(lon(aux),lat(aux),lags,1/1200,'dtime',60);
%     end
% end

vs(ni,:)=v(:,1)';
phis(ni,:)=v(:,2)';

%

end
t=toc
vvv{np}=vs;
phipp{np}=phis;
end
tt=(up+down)/2;
save([folder,'/prop_vel_corr_',data,'.mat'],'estacions','lon','lat','periods','limP','vvv','phipp','up','down','tt','corr_cell','mlag_cell','lw','DT')


%%
figure('Position',[50 50 1000 600])
colors=jet(length(periods));
pp=[]; pv=[];
for np=1:length(periods)
vs=vvv{np};
phis=phipp{np};

lm=0.09;
bm=0.08;
sp=0.001;
wid=0.88;
tp=0.04;
hi=(1-2*sp-bm-tp)/2;

if np==1
n=2;
loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi];
ax1=subplot('Position',loc);
hold on
tt=(up+down)/2;
else
end
% plot(ax1,tt,vs,'color',[0.8 0.8 0.8])
% hold on
% p1=plot(tt,prctile(vs,25,2),'-.r','LineWidth',1);
% hold on
% p2=plot(tt,prctile(vs,75,2),'--r','LineWidth',1);
% hold on
pv=[pv plot(ax1,tt,nanmedian(vs,2),'color',colors(np,:),'LineWidth',2)];
aux=not(isnan(nanmedian(vs,2)));
patch(ax1,[tt(aux) flip(tt(aux))],[prctile(vs(aux,:),10,2);flip(prctile(vs(aux,:),90,2))],colors(np,:),'FaceAlpha',0.1,'LineStyle','none')

ylim([10 50])
ylabel('Propagation speed [m/s]')
grid on
xlim([tt(1) tt(end)])
% title(['Periode: ',num2str(periods(np))])

%xlim([datenum(2021,5,23,7,0,0) datenum(2021,5,24,3,0,0)])
%xlim([datenum(2021,6,18,20,0,0) datenum(2021,6,20,7,0,0)])
ax=gca;
xl=xlim;
ticks=xl(1):2/24:xl(end);
ax.XTick=ticks;
cs=cos(phis); sn=sin(phis);
mphi=atan(sum(sn,2,'omitnan')./sum(cs,2,'omitnan'));

phis_de=270-phis*180/pi; mphi=270-mphi*180/pi;

aux=phis_de<=0;
phis_de(aux)=360+phis_de(aux);

if np==1
n=1;
loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi];
ax2=subplot('Position',loc);
hold on
vs(:,7)=NaN*ones(length(up),1);
tt=(up+down)/2;
else 
end
% p0=plot(tt,phis_de,'color',[0.5 0.5 0.5]);
% hold on
% p1=plot(tt,prctile(phis_de,25,2),'-.','LineWidth',1,'color','b');
% hold on
% p2=plot(tt,prctile(phis_de,75,2),'--','LineWidth',1,'color','b');
% hold on
pp=[pp plot(tt,nanmedian(phis_de,2),'color',colors(np,:),'LineWidth',2)];

aux=not(isnan(nanmedian(phis_de,2)));
patch([tt(aux) flip(tt(aux))],[prctile(phis_de(aux,:),10,2);flip(prctile(phis_de(aux,:),90,2))],colors(np,:),'FaceAlpha',0.1,'LineStyle','none')
% p4=plot(tt,ones(length(tt),length(c)).*phi,'--g','LineWidth',1.5);
% for n=1:length(c)
% p4=plot(tt,ones(length(tt),1).*phi(n),'--','LineWidth',1.5,'color',colors(n,:));
% end
% p4=plot(time(auxt),alph_p,'-s','LineWidth',1,'MarkerEdgeColor','k','color','g')

% datetick


ylabel('Propagation direction')
grid on
ylim([130 300])
xlim([tt(1) tt(end)])
%xlim([datenum(2021,5,23,7,0,0) datenum(2021,5,24,3,0,0)])
%xlim([datenum(2021,6,18,20,0,0) datenum(2021,6,20,7,0,0)])
ax=gca;
xl=xlim;
ticks=xl(1):2/24:xl(end);
ax.XTick=ticks;
ax.XTickLabel=datestr(ticks','HH:MM');

% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end



end
per=cellstr(num2str(periods'));
legend(ax1,pv,per,'Location','Northeast')
legend(ax2,pp,per,'Location','southwest')

end

return
%--------------------------------------------------------------------------
%%

% la0=[39.2];
% lo0=[2.7];
% f=[1/2000]; % secs
% c=[35];
% gamma=[pi/4];
% sigma=[400000]; % metres
% width=[20000]; % metres
% delay=[100]*1000;
% 
% la0=[37 37 39.2 37];
% lo0=[0 0.5 2.7 3];
% f=[1/1200 1/3600 1/1400 1/2000]; % secs
% c=[40 30 22 35];
% gamma=[pi/4 pi/3 pi/4 pi/7];
% sigma=[200000 200000 100000 400000]; % metres
% width=[30000 50000 25000 50000]; % metres
% delay=[0 150 320 100]*1000;

lats=[38.0:0.04:40.5]; lons=[1:0.04:5];
map=cell(length(time),1);
figure
m_proj('mercator','long',[lons(1) lons(end)],'lat',[lats(1) lats(end)]);
hold on

for nt=1:length(tt)
   mat=zeros(length(lats),length(lons));
    for nx=1:length(lons)
        for ny=1:length(lats)
            mat(ny,nx)=ona(lons(nx),lats(ny),tt(nt));
        end
    end
    map{nt}=mat;
m_pcolor(lons,lats,mat)
m_gshhs_i('line','LineWidth',1.5);
m_grid('linestyle','none','linewidth',2,'tickdir','in','xaxisloc','top','fontsize',15,'yticklabel','','xticklabel','');
m_quiver(lon,lat,vs(nt,:).*cos(phis(nt,:)),vs(nt,:).*sin(phis(nt,:)),'color','k','LineWidth',1)
% saveas(gcf,['./gif_vel/map_',num2str(nt)],'png')
pause(0.01)
end


%%

figure
plot(time,P_prova)

[P,ff,th,DOF]=peridiograms_v1(time/60,P_prova,1,'pwelch','l_win',1024);

figure
semilogy(1./ff,P)
xlim([2 80])
grid on
%%
function [patm]=ona(lon,lat,t)
global gamma f RT lon0 lat0 c sigma delay width la0 lo0
y=(lat-lat0)*RT*pi/180;
x=(lon-lon0)*RT*pi/180*cos(lat*pi/180);

y0=(lat-la0)*RT*pi/180;
x0=(lon-lo0)*RT*pi/180*cos(lat*pi/180);
for nc=1:length(c)
beta=2*pi*f(nc)/c(nc);
k=beta*cos(gamma(nc));
l=beta*sin(gamma(nc));
if nc==1
%patm=exp(-((x-c(nc)*t*cos(gamma(nc))).^2+(y-c(nc)*t*sin(gamma(nc))).^2)./(sigma(nc)^2))%.*cos(k*x+l*y-t*f(nc)*2*pi);
dd2=(sin(gamma(nc)-atan(y0(nc)/x0(nc))))^2*(x0(nc).^2+y0(nc).^2);
%;
patm=exp(-dd2./(width(nc)^2)).*exp(-(((x-c(nc)*cos(gamma(nc))*t+delay(nc))/sin(gamma(nc))+(y-c(nc)*sin(gamma(nc))*t+delay(nc))/cos(gamma(nc))).^2)./(sigma(nc)^2)).*cos(k*x+l*y-t*f(nc)*2*pi)+0.1*randn(length(t),1);
% figure
% plot(patm)
else
pp=exp(-(((x-c(nc)*cos(gamma(nc))*t+delay(nc))/sin(gamma(nc))+(y-c(nc)*sin(gamma(nc))*t+delay(nc))/cos(gamma(nc))).^2)./(sigma(nc)^2)).*cos(k*x+l*y-t*f(nc)*2*pi);


dd2=(sin(gamma(nc)-atan(y0(nc)/x0(nc))))^2*(y0(nc).^2+x0(nc).^2);

patm=patm+pp.*exp(-dd2./(width(nc)^2));
% figure
% plot(pp)
end

end
end



