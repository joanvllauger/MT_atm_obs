addpath '../aux_functions/'
addpath '../m_map/'

%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
% close all
global gamma f RT lon0 lat0 c sigma delay width la0 lo0 noise mmm

noise=0.2;

factor=0.2;

gamma=pi/4; % rad
RT=6371000; % m

lat0=[37];
lon0=[0];
radi=20;
sig=0.1;

f=1/600; % secs
c=40;
sigma=10000; % metres
delay=0;


la0=[37.5 37 39.2 39.7 39.5];
lo0=[0 0.5 2.7 2.5 2.5];
f=[1/1200 1/3600 1/1400 1/2000 1/2500]; % secs
c=[40 30 22 35 27];
gamma=[pi/5 pi/3 pi/4 pi/7 pi/2.5];
sigma=[200000 200000 100000 400000 200000]; % metres
width=[30000 50000 25000 50000 40000]; % metres
delay=[0 150 320 100 200]*1000;


limP=1./(0.8*[1:2:10])*60;
limP=flip(2*2.^([0:5]*7/5));
limP=limP(2:end-1);
limP=[60 30 15 8];
mm=movmedian(limP,2);

periods=mm(2:end)+1*randn(1,length(mm)-1)
la0=[37.5]*ones(1,length(periods));
lo0=[0]*ones(1,length(periods));
f=1./(periods*60); % 1/secs
c=linspace(20,40,length(periods));
c=[25 30 35].*ones(1,length(periods));
gamma=[pi/5 pi/4 pi/3].*ones(1,length(periods));
sigma=[3 2 1]*[100000].*ones(1,length(periods)); % metres
width=[800000]*ones(1,length(periods)); % metres
delay=[100 200 450]*1000.*ones(1,length(periods));

% periods=[10 20 50];
% la0=[37.5]*ones(1,length(periods));
% lo0=[0]*ones(1,length(periods));
% f=1./(periods*60); % 1/secs
% c=[40 25 30];
% gamma=[pi/5]*ones(1,length(periods));
% sigma=[100000 200000 900000];%*ones(1,length(periods)); % metres
% width=[80000]*ones(1,length(periods)); % metres
% delay=[300]*1000*2*ones(1,length(periods));



phi=270-gamma/pi*180
time=[0:60:43200]';
%--- Importam les posicions
load('../data/Atm_pres_all','mareografs','lon','lat')
estacions=mareografs;


mmm=0;
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
mmm=sqrt(mean(stds));

%--- Cream Patm
Patm=zeros(length(time),length(estacions));
for ne=1:length(estacions)
    Patm(:,ne)=ona(lon(ne),lat(ne),time);
end

for ne=[1 76]%:length(estacions)
figure
plot(time,Patm(:,ne))
title(estacions{ne})

figure
[P,ff,conf,DOF]=welch_v2(Patm(:,ne),512,1);
loglog(1./ff,P);
ax=gca;
ax.XTick=flip(periods);
ax.XTickLabel=num2str(flip(periods)');
grid on
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

df=3;
P=Patm;

vvv=cell(length(periods)+1,1);
phipp=cell(length(periods)+1,1);

stds=zeros(length(periods)+1,length(estacions));

for np=1:length(periods)+1

    if np<=length(periods)
    for ne=1:length(estacions)
    P(:,ne)=F2_filt_simple(time/60/60/24,Patm(:,ne),limP(np+1),limP(np));
    end
    else
        P=Patm;
    end
    
    %--- Eliminem les estacions que no no arriben acert valor de la
    % variança. 
    stds(np,:)=std(P,1,'omitnan');
    
    
    
    nn=30;
    figure
    
    plot(time,P(:,nn))
    hold on
    plot(time,Patm(:,nn))
    title(estacions{nn})
    
    P_prova=Patm(:,nn);
for ni=1:length(up)
    aux=time>=down(ni) & time<up(ni);

    P_t=P(aux,:);
    
    %--- VAriança d'aquest tram
    sd=std(P_t,1,'omitnan');
    aux=stds(np,:)*factor>sd;
    
    P_t(:,aux)=NaN;
    
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
         v(ne,:)=orlic_v_v1(lon(aux),lat(aux),est_corr_max(aux,aux),est_corr_mlag(aux,aux),sig,'dtime',60,'corr_th',0.9);
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

tt=down+(up-down)./2;
%% Dibuixam

ne=1;

N=4;

lm=0.07;
bm=0.06;
sp=0.0;
wid=0.86;
tp=0.05;
sp=0.04;

hi2=(1-bm-tp-sp*(N-1))/N;

hi=hi2/2;
cc=0.15/N;

figure('units','centimeters','Position',[2 2 17 15])
nn=1;
lletres=flip({'b','c','d'});
yl=3
for np=1:4
if np<4
    P=F2_filt_simple(time/60/60/24,Patm(:,ne),limP(np+1),limP(np));
else
    P=Patm(:,ne);
end

%--------------------------------------------------------------------------
n=2;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi)+cc,wid,hi-cc];
loc=[lm, bm+(sp+hi2)*(nn-1)+(n-1)*hi+cc,wid,hi-cc];
ax0=subplot('Position',loc);
plot(time,P)

% title(ax0,['(',lletres{nn},') ',datestr(times(nd),'dd-mmm-yyyy')])
ax0.TitleHorizontalAlignment='left';
% yl=max(abs(Patm(:,ne)));

% datetick
ax=gca;
% xl=xlim;
% ticks=xl(1):2/24:xl(end);
% ax.XTick=ticks;
ax.XTickLabel='';

grid on
xlim([time(1) time(end)])
if np<4
title(['(',lletres{np},') Filtered signal around period: ',num2str(periods(np),'%2.0f'),' minutes'])
ylim([-yl yl])
    ax.YTick=[-2:2:2];

else
    ylim([-5 5])
ax.YTick=[-4:4:4];
title('(a) Synthetic signal in Sa Rapita (South Mallorca)')
end
%---------------------------------------------------

pp=[]; pv=[];

% vvv=vvv(1:end-1);
% phis=phipp(1:end-1);

% np=1;

phis=phipp{np};
vs=vvv{np};



n=1;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi+cc];
loc=[lm, bm+(sp+hi2)*(nn-1)+(n-1)*hi,wid,hi+cc];
ax1=subplot('Position',loc);
hold on
if np<4
p4=plot(time,ones(length(time),1).*c(np),'--','LineWidth',1.5,'color','b');
else
  p4=plot(time,ones(length(time),1).*c,'--','LineWidth',1,'color','b');
  
end
pv=[pv plot(ax1,tt,nanmedian(vs,2),'b','LineWidth',2)];


aux=not(isnan(nanmedian(vs,2)));
patch(ax1,[tt(aux) flip(tt(aux))],[prctile(vs(aux,:),10,2);flip(prctile(vs(aux,:),90,2))],'b','FaceAlpha',0.1,'LineStyle','none')

ylim([15 45])
ylabel('Prop. speed [m/s]')
grid on
xlim([time(1) time(end)])
ax1.YTick=[20:5:40];
% title(ax1,datestr(times(nd),'dd-mmm-yyyy'))

%xlim([datenum(2021,5,23,7,0,0) datenum(2021,5,24,3,0,0)])
% xlim([times(nd)-wint/2 times(nd)+wint/2])
% xlim([t_in(nd) t_out(nd)])

% xl=xlim;
% ticks=xl(1):2/24:xl(end);
% ax1.XTick=ticks;
ax1.XTickLabel='';
ax1.Visible='on';


cs=cos(phis); sn=sin(phis);
mphi=atan(sum(sn,2,'omitnan')./sum(cs,2,'omitnan'));

phis_de=270-phis*180/pi; mphi=270-mphi*180/pi;

aux=phis_de<=0;
phis_de(aux)=360+phis_de(aux);

% n=1;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi];
ax2=axes('Position',loc);
hold on
vs(:,7)=NaN*ones(length(tt),1);

% p0=plot(tt,phis_de,'color',[0.5 0.5 0.5]);
% hold on
% p1=plot(tt,prctile(phis_de,25,2),'-.','LineWidth',1,'color','b');
% hold on
% p2=plot(tt,prctile(phis_de,75,2),'--','LineWidth',1,'color','b');
% hold on
if np<4
p4=plot(time,ones(length(time),1).*phi(np),'--','LineWidth',1.5,'color','r');
else
    p4=plot(time,ones(length(time),1).*phi,'--','LineWidth',1,'color','r');

end
pp=[pp plot(ax2,tt,nanmedian(phis_de,2),'r','LineWidth',2)];
aux=not(isnan(nanmedian(phis_de,2)));
patch([tt(aux) flip(tt(aux))],[prctile(phis_de(aux,:),10,2);flip(prctile(phis_de(aux,:),90,2))],'r','FaceAlpha',0.1,'LineStyle','none')
% p4=plot(tt,ones(length(tt),length(c)).*phi,'--g','LineWidth',1.5);
% for n=1:length(c)
% p4=plot(tt,ones(length(tt),1).*phi(n),'--','LineWidth',1.5,'color',colors(n,:));
% end
% p4=plot(time(auxt),alph_p,'-s','LineWidth',1,'MarkerEdgeColor','k','color','g')

% datetick


ylabel('Prop. direction [º]')
ax2.YAxisLocation='right';
ax2.Color='none';
grid on
ylim([200 260])

% xlim([times(nd)-wint/2 times(nd)+wint/2])
% xlim([t_in(nd) t_out(nd)])
% xl=xlim;
% ticks=xl(1):2/24:xl(end);
% ax2.XTick=ticks;
% ax2.XTickLabel=datestr(ticks','HH:MM');
ax2.YTick=[210:10:250];
% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end
xlim([time(1) time(end)])


% 
% plot(ax1,ww.time,ww.wspeed,'k','LineWidth',1,'Marker','s','MarkerFaceColor','g')
% plot(ax2,ww.time,270-ww.wdir*180/pi,'k','LineWidth',1,'Marker','s','MarkerFaceColor','g')

% per=strcat(num2str(round(limP(2:end-1)')),' - ');
% per=strcat(per,num2str(round(limP(1:end-2)')));

% per=cellstr(strcat(per,' min'));
ax2.XTick=[0:12]*60*60;
ax0.XTick=[0:12]*60*60;
ax1.XTick=[0:12]*60*60;
% legend(ax1,pv,per,'Location','Northeast')
if nn==1
    ax2.XTickLabel=num2str(ax2.XTick'/60/60,'%2.0f');
    xlabel('Time [h]')
else
    ax2.XTickLabel='';
end

if nn==3
legend(ax2,[pv pp],{'Propagation speed','Propagation direction'},'Location','northwest')
end


nn=nn+1;

end

% ax4=axes('Position',[lm 1-tp+sp wid tp-2*sp]);
% plot(time,Patm(:,ne),'r')
% ax4.XTickLabel='';
% xlim([time(1) time(end)])
% title('(a) Synthetic signal in Sa Rapita (South Mallorca)')
% grid on
% ax4.TitleHorizontalAlignment='left';
% ylabel(['hPa'])
% ax4.XTick=[0:12]*60*60;
print(gcf,'-djpeg','-r500',['SI1_5'])
%%

% [vs,IA,IC]=unique(V','rows'); vs=vs';
% lon=LON(IA); lat=LAT(IA); phis=phis(:,IA); estacions=estacions(IA);
% aux=lon<2;
% 
% vs(:,aux)=[]; phis(:,aux)=[];
% lon(aux)=[]; lat(aux)=[]; estacions(aux)=[];
% vs=V;
% phis=PHI;
for np=1:length(periods)
vs=vvv{np};
phis=phipp{np};

lm=0.09;
bm=0.08;
sp=0.001;
wid=0.88;
tp=0.04;
hi=(1-2*sp-bm-tp)/2;

figure('Position',[50 50 1000 600])
n=2;
loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi];
subplot('Position',loc)
tt=(up+down)/2;
plot(tt,vs,'color',[0.8 0.8 0.8])
hold on
p1=plot(tt,prctile(vs,25,2),'-.r','LineWidth',1);
hold on
p2=plot(tt,prctile(vs,75,2),'--r','LineWidth',1);
hold on
p3=plot(tt,nanmedian(vs,2),'k','LineWidth',2);

colors=jet(length(c));
for n=1:length(c)
p4=plot(tt,ones(length(tt),1).*c(n),'--','LineWidth',1.5,'color',colors(n,:));
end
% p4=plot(time(auxt),c_p,'-o','LineWidth',1,'MarkerEdgeColor','k','color','g');
ylim([10 50])
ylabel('Propagation speed [m/s]')
grid on
xlim([tt(1) tt(end)])
title(['Periode: ',num2str(periods(np))])
legend([p1 p2 p3],{'Percentile 25%','Percentile 75%','Median'},'Location','Northeast')

cs=cos(phis); sn=sin(phis);
mphi=atan(sum(sn,2,'omitnan')./sum(cs,2,'omitnan'));

phis_de=270-phis*180/pi; mphi=270-mphi*180/pi;

aux=phis_de<=0;
phis_de(aux)=360+phis_de(aux);

n=1;
loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi];
subplot('Position',loc)
vs(:,7)=NaN*ones(length(up),1);
tt=(up+down)/2;
p0=plot(tt,phis_de,'color',[0.5 0.5 0.5]);
hold on
p1=plot(tt,prctile(phis_de,25,2),'-.','LineWidth',1,'color','b');
hold on
p2=plot(tt,prctile(phis_de,75,2),'--','LineWidth',1,'color','b');
hold on
p3=plot(tt,nanmedian(phis_de,2),'k','LineWidth',2);
% p4=plot(tt,ones(length(tt),length(c)).*phi,'--g','LineWidth',1.5);
for n=1:length(c)
p4=plot(tt,ones(length(tt),1).*phi(n),'--','LineWidth',1.5,'color',colors(n,:));
end
% p4=plot(time(auxt),alph_p,'-s','LineWidth',1,'MarkerEdgeColor','k','color','g')




ylabel('Propagation direction')
grid on
ylim([180 270])
xlim([tt(1) tt(end)])
% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end
legend([p1 p2 p3],{'Percentile 25%','Percentile 75%','Median'},'Location','southwest')


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
global gamma f RT lon0 lat0 c sigma delay width la0 lo0 noise mmm 
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
patm=exp(-dd2./(width(nc)^2)).*exp(-(((x-c(nc)*cos(gamma(nc))*t+delay(nc))/sin(gamma(nc))+(y-c(nc)*sin(gamma(nc))*t+delay(nc))/cos(gamma(nc))).^2)./(sigma(nc)^2)).*cos(k*x+l*y-t*f(nc)*2*pi);
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

mmm=std(patm);
patm=patm+mmm*noise*randn(length(t),1);
end






