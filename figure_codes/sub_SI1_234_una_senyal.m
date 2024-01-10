
%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
% close all
global gamma f RT lon0 lat0 c sigma delay width la0 lo0 noise dt phase mm

RT=6371000; % m

lat0=[37];
lon0=[0];
radi=20;
sig=0.1;

la0=[37.5]*ones(1,length(periods));
lo0=[0]*ones(1,length(periods));
f=1./(periods*60); % 1/secs
c=vel*ones(1,length(periods));
gamma=deg*ones(1,length(periods));
sigma=[1000000]*ones(1,length(periods)); % metres
width=[800000]*ones(1,length(periods)); % metres
delay=[150]*1000*ones(1,length(periods));
phase=rand(length(periods),1);

phi=270-gamma/pi*180;
time=[0:dt:43200]';
%--- Importam les posicions
load('../data/Atm_pres_all','mareografs','lon','lat')
estacions=mareografs;

mm=0;
lats=[38.0:0.1:40.5]; lons=[1:0.1:5];
stds=zeros(length(time),1);
for nt=1:length(time)
   mat=zeros(length(lats),length(lons));
    for nx=1:length(lons)
        for ny=1:length(lats)
            mat(ny,nx)=ona(lons(nx),lats(ny),time(nt));
        end
    end
    stds(nt)=var(mat,1,'all');

end
mm=sqrt(mean(stds));


%--- Cream Patm
Patm=zeros(length(time),length(estacions));
for ne=1:length(estacions)
    Patm(:,ne)=ona(lon(ne),lat(ne),time);
end

if k==3 && n==4
for ne=[1 76]%:length(estacions)
figure
plot(time,Patm(:,ne))
title(estacions{ne})

figure
[P,ff,conf,DOF]=welch_v2(Patm(:,ne),512,1);
loglog(1./ff,P);
ax=gca;
ax.XTick=periods;
ax.XTickLabel=num2str(periods')
grid on
end
end


%%
%--- time intervals
tini=time(1); tfi=time(end);
DT=30*60;
lw=2*60*60;
up=[tini+lw:DT:tfi];
down=[tini:DT:tfi-lw];



%--- arreglam le sposicions
aux=isnan(lon)|isnan(lat); 
lat(aux)=[]; lon(aux)=[]; Patm(:,aux)=[]; estacions(aux)=[];

%-------------------------------------------------------------------------- 
%%
%--- definim els vetor per guardar les velocitats
vs=ones(length(up),length(estacions));
phis=vs;
std_mat=ones(length(up),1);

%--- Mesurem lesdistàncies
Ne=length(lat);
Dist=zeros(Ne,Ne);
for ne=1:Ne
Dist(ne,:)=deg2km(distance(lat(ne),lon(ne),lat,lon));
end

%%
freq=10;
df=7.5;
%--- Col·locació del figures
lm=0.07;
bm=0.07;
sp=0.07;
wid=(1-2*lm-0.03)/2;
hi=(1-3*sp-bm)/3;
%%
tic
for ni=1:length(up)
    aux=time>=down(ni) & time<=up(ni);

    P_t=Patm(aux,:);
%     for ne=1:length(estacions)
%     P_t(:,ne)=F2_filt_simple(time,P_t(:,ne),freq-df,freq+df);
%     end

    

[est_corr_max,est_corr_mlag]=F3_correlation_v6(P_t,Dist,radi);


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
         v(ne,:)=orlic_v_v1(lon(aux),lat(aux),est_corr_max(aux,aux),est_corr_mlag(aux,aux),sig,'dtime',60);
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

std_mat(ni)=mean(nanstd(P_t,[],1));
%

end
t=toc


%% Dibuixam
cs=cos(phis); sn=sin(phis);
mphi=atan(sum(sn,2,'omitnan')./sum(cs,2,'omitnan'));
phis_de=270-phis*180/pi; mphi=270-mphi*180/pi;

aux=phis_de<=0;
phis_de(aux)=360+phis_de(aux);

if k==3 && n==4
% [vs,IA,IC]=unique(V','rows'); vs=vs';
% lon=LON(IA); lat=LAT(IA); phis=phis(:,IA); estacions=estacions(IA);
% aux=lon<2;
% 
% vs(:,aux)=[]; phis(:,aux)=[];
% lon(aux)=[]; lat(aux)=[]; estacions(aux)=[];
% vs=V;
% phis=PHI;

lm=0.09;
bm=0.08;
sp=0.001;
wid=0.88;
tp=0.28;
hi=(1-2*sp-bm-tp)/2;

figure('units','centimeters','Position',[2 2 17 10])
ne=6;
loc=[lm,1-tp,wid,0.25];
ax3=axes('Position',loc);
plot(time/60/60,Patm(:,ne))
xlim([0 12])
grid on
ax3.XTickLabel='';
% title(estacions{ne})
ylabel('Atm. Pres. [hPa]')
annotation('textbox',[lm loc(2)+0.25-0.06 0.5 0.06],'String','(a) Synthetic signal (Ciutadella)','BackgroundColor','none','LineStyle','none','FontWeight','bold')


nnn=2;
loc=[lm,(bm+sp*(nnn-1)+(nnn-1)*hi),wid,hi];
subplot('Position',loc)
tt=(up+down)/2;
plot(tt/60/60,vs,'color',[0.8 0.8 0.8])
hold on
p1=plot(tt/60/60,prctile(vs,25,2),'-.r','LineWidth',1);
hold on
p2=plot(tt/60/60,prctile(vs,75,2),'--r','LineWidth',1);
hold on
p3=plot(tt/60/60,nanmedian(vs,2),'k','LineWidth',2);

colors=jet(length(c));
for nnn=1:length(c)
p4=plot([0 12],ones(2,1).*c(nnn),'-','LineWidth',1.2,'color','b');
end
% p4=plot(time(auxt),c_p,'-o','LineWidth',1,'MarkerEdgeColor','k','color','g');
ylim([25 35])
annotation('textbox',[lm loc(2)+hi-0.06 0.3 0.06],'String','(b) Propagation speed','BackgroundColor','none','LineStyle','none','FontWeight','bold')

ylabel('m/s')
grid on
xlim([0 12])
ax1=gca;
ax1.YTick=[26:1:34];
ax1.XTickLabel='';

% legend([p1 p2 p3],{'Percentile 25%','Percentile 75%','Median'},'Location','Northeast')





nnn=1;
loc=[lm,(bm+sp*(nnn-1)+(nnn-1)*hi),wid,hi];
subplot('Position',loc)
vs(:,7)=NaN*ones(length(up),1);
tt=(up+down)/2;
p0=plot(tt/60/60,phis_de,'color',[0.5 0.5 0.5]);
hold on
p1=plot(tt/60/60,prctile(phis_de,25,2),'-.','LineWidth',1,'color','r');
hold on
p2=plot(tt/60/60,prctile(phis_de,75,2),'--','LineWidth',1,'color','r');
hold on
p3=plot(tt/60/60,nanmedian(phis_de,2),'k','LineWidth',2);
% p4=plot(tt,ones(length(tt),length(c)).*phi,'--g','LineWidth',1.5);
for nnn=1:length(c)
p4=plot([0 12],ones(2,1).*phi(nnn),'-','LineWidth',1.2,'color','b');
end
% p4=plot(time(auxt),alph_p,'-s','LineWidth',1,'MarkerEdgeColor','k','color','g')

annotation('textbox',[lm loc(2)+hi-0.06 0.3 0.06],'String','(c) Propagation direction','BackgroundColor','none','LineStyle','none','FontWeight','bold')


ylabel('Degrees')
grid on
ylim([235 245])
xlim([0 12])
ax2=gca;
ax2.YTick=[235:1:244];
xlabel('Hours')
% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end
legend1=legend([p1 p2 p3 p4],{'Percentile 25%','Percentile 75%','Median','Theoretical value'},'Location','southwest')
set(legend1,...
    'Position',[0.754408405106591 0.38291161877804 0.208722737580074 0.135761585719801]);
print(gcf,'-djpeg','-r500',['SI1_2'])

end

%%
%-- Filtram
stdv=nanstd(vs,[],2);
menv=nanmean(vs,2);

stdp=nanstd(phis_de,[],2);
menp=nanmean(phis_de,2);

aux=vs>menv+3*stdv | vs<menv-3*stdv | phis_de>menp+3*stdp | phis_de<menp-3*stdp;
ver=vs;
per=phis_de;
% ver(aux)=NaN;
% per(aux)=NaN;

err=sqrt(nanmean((ver-vel).^2,1));
bias=nanmean((ver-vel),1);
errp=sqrt(nanmean((per-phi(1)).^2,1));
biasp=nanmean((per-phi(1)),1);

if 0
figure
plot(tt,std_mat,'o')
hold on
plot(tt,err,'o')
hold on
plot(tt,bias,'o')
title('v')
ylim([-3 3])
grid on
figure

plot(std_mat,err,'o')
ylabel('error')


figure
plot(tt,std_mat,'o')
hold on
plot(tt,errp,'o')
hold on
plot(tt,biasp,'o')
title('phi')
ylim([-3 3])

figure

plot(std_mat,errp,'o')
ylabel('error')
ylim([-3 3])
end

%%
if k==3 && n==4
figure('Units','centimeters','Position',[2 2 18 11.2])

MS=15;

lm=0.01;
bm=0.005;
vvs=0.02;
hs=0.1;
tp=0.03;
rm=0.1;
hi=(1-vvs-bm-tp)/2;
wi=(1-hs-lm-rm)/2;

N=16;
cmap=pmkmp(N,'TwoColors');

%---
lats=[38.4:0.04:40.2]; lons=[1:0.04:4.4];
nnn=1; 
kkk=1;
pos=[lm+(wi+hs)*(kkk-1) 1-tp-hi-(hi+vvs)*(nnn-1) wi hi];
ax1=axes('Position',pos);

m_proj('mercator','long',[lons(1) lons(end)],'lat',[lats(1) lats(end)]);
m_gshhs_i('line','LineWidth',1,'color','k');
hold on
aux=not(isnan(err));
m_scatter(lon(aux),lat(aux),MS,err(aux),'filled','MarkerEdgeColor','k')
cb=colorbar;
cb.Label.String='m/s';
cb.Label.FontSize=10;
cb.Position(3)=cb.Position(3)*0.7;
cb.Position(1)=pos(1)+pos(3)+0.01;
cb.Position(4)=hi-0.04;
cb.Position(2)=pos(2)+0.02;
colormap(ax1,cmap(end-N/2+1:end,:))
caxis([0 2])
% colormap('jet')

m_grid('linestyle','none','linewidth',2,'tickdir','in','xaxisloc','top','fontsize',15,'yticklabel','','xticklabel','');
title('(a) RMSE Speed estimates')
% saveas(gcf,['./gif_vel/map_',num2str(nt)],'png')

%---
nnn=2; 
kkk=1;
pos=[lm+(wi+hs)*(kkk-1) 1-tp-hi-(hi+vvs)*(nnn-1) wi hi];
ax2=axes('Position',pos);
m_proj('mercator','long',[lons(1) lons(end)],'lat',[lats(1) lats(end)]);
m_gshhs_i('line','LineWidth',1,'color','k');
hold on
aux=not(isnan(errp));
m_scatter(lon(aux),lat(aux),MS,errp(aux),'filled','MarkerEdgeColor','k')
cb=colorbar;
cb.Label.String='Degrees (º)';
cb.Label.FontSize=10;
cb.Position(3)=cb.Position(3)*0.7;
cb.Position(1)=pos(1)+pos(3)+0.01;
cb.Position(4)=hi-0.04;
cb.Position(2)=pos(2)+0.02;
caxis([0 4])
colormap(ax2,cmap(end-N/2+1:end,:))
% colormap('jet')

m_grid('linestyle','none','linewidth',2,'tickdir','in','xaxisloc','top','fontsize',15,'yticklabel','','xticklabel','');
title('(c) RMSE Direction estimates')
% saveas(gcf,['./gif_vel/map_',num2str(nt)],'png')


nnn=1; 
kkk=2;
pos=[lm+(wi+hs)*(kkk-1) 1-tp-hi-(hi+vvs)*(nnn-1) wi hi];
ax3=axes('Position',pos);

m_proj('mercator','long',[lons(1) lons(end)],'lat',[lats(1) lats(end)]);
m_gshhs_i('line','LineWidth',1,'color','k');
hold on
aux=not(isnan(bias));
m_scatter(lon(aux),lat(aux),MS,bias(aux),'filled','MarkerEdgeColor','k')
cb=colorbar;
cb.Label.String='m/s';
cb.Label.FontSize=10;
cb.Position(3)=cb.Position(3)*0.7;
cb.Position(1)=pos(1)+pos(3)+0.01;
cb.Position(4)=hi-0.04;
cb.Position(2)=pos(2)+0.02;
caxis([-2 2])
colormap(ax3,cmap)
% colormap('jet')

m_grid('linestyle','none','linewidth',2,'tickdir','in','xaxisloc','top','fontsize',15,'yticklabel','','xticklabel','');
title('(b) BIAS Speed estimates')


nnn=2; 
kkk=2;
pos=[lm+(wi+hs)*(kkk-1) 1-tp-hi-(hi+vvs)*(nnn-1) wi hi];
ax4=axes('Position',pos);

m_proj('mercator','long',[lons(1) lons(end)],'lat',[lats(1) lats(end)]);
m_gshhs_i('line','LineWidth',1,'color','k');
hold on
aux=not(isnan(biasp));
m_scatter(lon(aux),lat(aux),MS,biasp(aux),'filled','MarkerEdgeColor','k')
cb=colorbar;
cb.Label.String='Degrees (º)';
cb.Label.FontSize=10;
cb.Position(3)=cb.Position(3)*0.7;
cb.Position(1)=pos(1)+pos(3)+0.01;
cb.Position(4)=hi-0.04;
cb.Position(2)=pos(2)+0.02;
caxis([-4 4])
colormap(ax4,cmap)

m_grid('linestyle','none','linewidth',2,'tickdir','in','xaxisloc','top','fontsize',15,'yticklabel','','xticklabel','');
title('(d) BIAS Direction estimates')

ax1.TitleHorizontalAlignment = 'left';
ax2.TitleHorizontalAlignment = 'left';
ax3.TitleHorizontalAlignment = 'left';
ax4.TitleHorizontalAlignment = 'left';
% 
print(gcf,'-djpeg','-r500',['SI1_3'])
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
% [patm]=ona(1,38,time);
function [patm]=ona(lon,lat,t)
global gamma f RT lon0 lat0 c sigma delay width la0 lo0 noise dt phase mm
y=(lat-lat0)*RT*pi/180;
x=(lon-lon0)*RT*pi/180*cos(lat*pi/180);

y0=(lat-la0)*RT*pi/180;
x0=(lon-lo0)*RT*pi/180*cos(lat*pi/180);
for nc=1:length(c)
beta=2*pi*f(nc)/c(nc)*dt/60;
k=beta*cos(gamma(nc));
l=beta*sin(gamma(nc));
if nc==1
%patm=exp(-((x-c(nc)*t*cos(gamma(nc))).^2+(y-c(nc)*t*sin(gamma(nc))).^2)./(sigma(nc)^2))%.*cos(k*x+l*y-t*f(nc)*2*pi);
%dd2=(sin(gamma(nc)-atan(y0(nc)/x0(nc))))^2*(x0(nc).^2+y0(nc).^2);
%;
%patm=exp(-dd2./(width(nc)^2)).*
%exp(-(((x-c(nc)*cos(gamma(nc))*t+delay(nc))/sin(gamma(nc))+(y-c(nc)*sin(gamma(nc))*t+delay(nc))/cos(gamma(nc))).^2)./(sigma(nc)^2));
patm=cos(k*x+l*y-t*f(nc)*2*pi+2*pi*phase(nc));
else
%exp(-(((x-c(nc)*cos(gamma(nc))*t+delay(nc))/sin(gamma(nc))+(y-c(nc)*sin(gamma(nc))*t+delay(nc))/cos(gamma(nc))).^2)./(sigma(nc)^2))
pp=cos(k*x+l*y-t*f(nc)*2*pi+2*pi*phase(nc));


%dd2=(sin(gamma(nc)-atan(y0(nc)/x0(nc))))^2*(y0(nc).^2+x0(nc).^2);

patm=patm+pp;% .*exp(-dd2./(width(nc)^2));
% figure
% plot(pp)
end

end

patm=patm+mm*noise*randn(length(t),1);
end