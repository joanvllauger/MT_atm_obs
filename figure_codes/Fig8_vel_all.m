addpath '../aux_functions/'
% addpath C:\Users\Usuario\Desktop\VENOM\Codis\Analisi\Article_atm\prop_vel

%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
folder=dir('..\meteotsunamis\meteo*');
files={folder(:).name}';
dates=cellfun(@(x) strrep(x,'meteotsunami_',''),files,'UniformOutput',0);
wint=12/24; % dies
add={'','','','','','a','b','','',''};


sl=load('../data/ciutadella_SL_AtmPres.mat','time','SL');
%%
%--- Importam les rissagues que tenim controlades
T=readtable('../list_of_events.xlsx');
ww=load('../data/wind_ciutadella.mat');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));
wh=datenum(table2array(T(:,4)));
%%
lm=0.09;
bm=0.08;
sp=0.01;
wid=0.88;
tp=0.04;
hi=(1-3*sp-bm-tp)/3;

MS=70;
MS=7;
fact=25;
DT=6/24;

dd=75; 
uu=25;



x=zeros(length(times),5);
y=x;
yn=x;
yp=x;
xn=x;
xp=x;

for nd=1:length(times)
    datestr(times(nd))
    [mm,aaa]=min(abs(datenum(dates,'yyyymmdd')-times(nd)));
    data=dates{aaa};
    data=datestr(times(nd),'yyyymmdd');
  folder=['../meteotsunamis/meteotsunami_',data,'/'];
    
load([folder,'prop_vel_corr_',data,'.mat'],'estacions','periods','limP','vvv','phipp','tt')

pp=[]; pv=[];

vvv=vvv(1:end-1);
phis=phipp(1:end-1);
phis=cellfun(@(x) 270-x*180/pi,phis,'UniformOutput',0);
periods=periods(1:end-1);
aux=tt>=times(nd)-DT & tt<=times(nd);


% h1=errorbar(nanmedian(phis{2}(aux,:),'all'),nanmedian(vvv{2}(aux,:),'all'),prctile(vvv{2}(aux,:),dd,'all'),prctile(vvv{2}(aux,:),uu,'all'),prctile(phis{2}(aux,:),dd,'all'),prctile(phis{2}(aux,:),uu,'all'),'MarkerSize',MS,'MarkerFaceColor',colors(nd,:),'Marker','s','MarkerEdgeColor','k');
% h2=errorbar(nanmedian(phis{3}(aux,:),'all'),nanmedian(vvv{3}(aux,:),'all'),prctile(vvv{3}(aux,:),dd,'all'),prctile(vvv{3}(aux,:),uu,'all'),prctile(phis{3}(aux,:),dd,'all'),prctile(phis{3}(aux,:),uu,'all'),'MarkerSize',MS,'MarkerFaceColor',colors(nd,:),'Marker','o','MarkerEdgeColor','k');
h=[];
for n=2:4
    x(nd,n-1)=nanmedian(phis{n}(aux,:),'all');
    y(nd,n-1)=nanmedian(vvv{n}(aux,:),'all');
    yn(nd,n-1)=y(nd,n-1)-prctile(vvv{n}(aux,:),dd,'all');
    yp(nd,n-1)=prctile(vvv{n}(aux,:),uu,'all')-y(nd,n-1);
    xn(nd,n-1)=x(nd,n-1)-prctile(phis{n}(aux,:),dd,'all');
    xp(nd,n-1)=prctile(phis{n}(aux,:),uu,'all')-x(nd,n-1);
end

load([folder,'prop_vel_corr_allfreq_',data,'.mat'],'estacions','vvv','phipp','tt')
phis=phipp;
phis=270-phis*180/pi;

aux=tt>=times(nd)-DT & tt<=times(nd);

x(nd,4)=nanmedian(phis(aux,:),'all');
y(nd,4)=nanmedian(vvv(aux,:),'all');
yn(nd,4)=y(nd,4)-prctile(vvv(aux,:),dd,'all');
yp(nd,4)=prctile(vvv(aux,:),uu,'all')-y(nd,4);
xn(nd,4)=x(nd,4)-prctile(phis(aux,:),dd,'all');
xp(nd,4)=prctile(phis(aux,:),uu,'all')-x(nd,4);
%--------------------------------------------------------------------------

phis=ww.wdir;
phis=270-phis*180/pi;
vvv=ww.wspeed;
aux=ww.time>=times(nd)-DT & ww.time<=times(nd);

x(nd,5)=nanmedian(phis(aux),'all');
y(nd,5)=nanmedian(vvv(aux),'all');
yn(nd,5)=y(nd,5)-prctile(vvv(aux),dd,'all');
yp(nd,5)=prctile(vvv(aux),uu,'all')-y(nd,5);
xn(nd,5)=x(nd,5)-prctile(phis(aux),dd,'all');
xp(nd,5)=prctile(phis(aux),uu,'all')-x(nd,5);
end
%%
N=2;
lm=0.08;
rm=0.08;
tp=0.02;
sp=0.03;
bm=0.05;

wi=1-(lm+rm);
hi=(1-(tp+bm+sp*(N-1)))/2;


xl=[ 0.5 10.5];

MS=5;
ttt=1:length(times);
figure('units','centimeters','Position',[2 2 20 12])
markers={'o','o','o','d','v'};
colors=hsv(3);
colors=[colors; [1 1 1]*0; [1 1 0]*1];

ax1=subplot(211,'Position',[lm 1-tp-hi wi hi]);
hold on
patch([0 0 xl(2) xl(2)],[24 36 36 24]*1/50+0.3,[1 1 1]*0.9,'LineStyle','none')
plot(ttt,wh,'s','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS*1.4,'LineWidth',1.2)
ax1.XLim=xl;
% ax1.Visible = 'off';
ax1.YAxisLocation='right';
tics=[0.3:0.1:1.3]';
ax1.Visible='on';
ax1.XTickLabel='';
ax1.YTick=tics;
ylim([0.3 tics(end)])
ax1.YGrid='off';
ylabel(ax1,'Maximum SLO height [m]');
% ax1.GridColor='r';
% ax1.GridAlpha=0.5;
loc=ax1.Position;
xlim(xl)

ax1.YTickLabel=num2str(tics);

ax=axes('Position',loc);
hold on
h=[];
for n=2:6
    if n==5
        h=[h errorbar(ttt+(4.5-n)*0.15,y(:,n-1),yn(:,n-1),yp(:,n-1),'MarkerSize',MS*1.3,'Marker',markers{n-1},'MarkerEdgeColor','k','MarkerFaceColor',colors(n-1,:),'LineStyle','none','color','k')];
    elseif n==6
        h=[h errorbar(ttt+(4.5-n)*0.15,y(:,n-1),yn(:,n-1),yp(:,n-1),'MarkerSize',MS*1,'Marker',markers{n-1},'MarkerEdgeColor','k','MarkerFaceColor',colors(n-1,:),'LineStyle','none','color','k')];
 
    else
        h=[h errorbar(ttt+(4.5-n)*0.15,y(:,n-1),yn(:,n-1),yp(:,n-1),'MarkerSize',MS*0.7,'Marker',markers{n-1},'MarkerEdgeColor','k','MarkerFaceColor',colors(n-1,:),'LineStyle','none','color',colors(n-1,:))];
    end
end
ylim([0 50])
xlim([ttt(1)-1 ttt(end)+1])
grid on
ax=gca;
ylabel(ax,'Propagation speed [m/s]')

ax.YTick=[0:5:50];
% xl=ax.XLim;
ax.XTickLabel='';
ax.Color='none';
ax.XGrid='off';
leg=strcat(num2str(limP(3:end-1)','%3.0f'),' - ');
leg=strcat(leg,num2str(limP(2:end-2)','%3.0f'));
leg=strcat(leg,' min');
leg=cellstr(leg);
leg=cat(1,leg,{'2 - 120 min'});
leg=cat(1,leg,{'ERA5 500hPa'});
legend(flip(h),flip(leg),'Location','southwest','orientation','horizontal')
xlim(xl)

ax1=subplot(212,'Position',[lm bm wi hi])
patch([0 0 xl(2) xl(2)],([210 260 260 210]-150)*1/150+0.3,[1 1 1]*0.9,'LineStyle','none')
hold on
% patch([0 0 10 10],[24 36 36 24]*2/50-0.5,[1 1 1]*0.9,'LineStyle','none')
% plot(ttt,wh,'-v','color',[1 1 1]*0.5,'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]*0.5,'MarkerSize',MS)
ax1.XLim=xl;
% ax1.Visible = 'off';
ax1.YAxisLocation='right';
tics=[0.3:0.1:1.3]';
ax1.Visible='on';
ax1.XTickLabel='';
ax1.YTick=tics;
ax1.YTickLabel='';
ylim([0.3 tics(end)])
ax1.YGrid='off';
% ylabel(ax1,'Maximum SLO height [m]');
% ax1.GridColor='r';
% ax1.GridAlpha=0.5;
loc=ax1.Position;
xlim(xl)

ax1.Color='w';
ax1.XTick=ttt;

ssss=cat(2,cellstr(datestr(t_in,'mmm-dd')),cellstr(add)');
sssss=cellfun(@(aa,bb) [aa,bb],ssss(:,1),ssss(:,2),'UniformOutput',0);

ax1.XTickLabel=sssss;

loc=ax1.Position;

ax=axes('Position',loc);
hold on
h=[];
for n=2:6
    if n==5
        h=[h errorbar(ttt+(4.5-n)*0.1,x(:,n-1),xn(:,n-1),xp(:,n-1),'MarkerSize',MS*1.3,'Marker',markers{n-1},'MarkerEdgeColor','k','MarkerFaceColor',colors(n-1,:),'LineStyle','none','color',colors(n-1,:))];
    elseif n==6
        h=[h errorbar(ttt+(4.5-n)*0.1,x(:,n-1),xn(:,n-1),xp(:,n-1),'MarkerSize',MS,'Marker',markers{n-1},'MarkerEdgeColor','k','MarkerFaceColor',colors(n-1,:),'LineStyle','none','color','k')];
    else
h=[h errorbar(ttt+(4.5-n)*0.1,x(:,n-1),xn(:,n-1),xp(:,n-1),'MarkerSize',MS*0.7,'Marker',markers{n-1},'MarkerEdgeColor','k','MarkerFaceColor',colors(n-1,:),'LineStyle','none','color',colors(n-1,:))];
    end
end
xlim([ttt(1)-1 ttt(end)+1])
ylim([150 300])
grid on
ylabel('Propagation direction [º]')

ax=gca;
ax.XTickLabelRotation=40;
ax.XTickLabel='';
ax.YTick=[150:15:300]
loc=ax.Position;

ax.Color='none';
ax.XGrid='off';
xlim(xl)

t=annotation('textbox',[lm bm+hi-0.05 0.1 0.05],'String','(b)','VerticalAlignment','bottom');
t.FontSize=10;
t.FontWeight='bold';
t.LineStyle='none';

t=annotation('textbox',[lm bm+hi+sp+hi-0.05 0.1 0.05],'String','(a)','VerticalAlignment','bottom');
t.FontSize=10;
t.FontWeight='bold';
t.LineStyle='none';



print(gcf,'-djpeg','-r500',['Fig8'])
return
%%
ddd=phis{3}(aux,:)
figure
histogram(ddd)
hold on
plot(nanmedian(ddd,'all')+ylim*0,ylim)
plot(prctile(ddd,25,'all')+ylim*0,ylim)
plot(prctile(ddd,75,'all')+ylim*0,ylim)