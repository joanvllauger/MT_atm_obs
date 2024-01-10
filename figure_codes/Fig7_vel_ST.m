addpath '../aux_functions/'
%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
folder=dir('..\meteotsunamis\meteo*');
files={folder(:).name}';
dates=cellfun(@(x) strrep(x,'meteotsunami_',''),files,'UniformOutput',0);
wint=12/24; % dies

sl=load('../data/ciutadella_SL_AtmPres.mat','time','SL');
%%
%--- Importam les rissagues que tenim controlades
T=readtable('../List_of_events.xlsx');
% ww=load('wind_ciutadella.mat');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));

%%
N=3;

lm=0.07;
bm=0.04;
sp=0.0;
wid=0.86;
tp=0.03;
sp=0.08;

hi2=(1-bm-tp-sp*(N-1))/N;

hi=hi2/2;
cc=0.15/N;

figure('units','centimeters','Position',[2 2 17 16])
nn=1;
lletres=flip({'a','b','c'});
%%
for nd=flip([3 4 9])
    datestr(times(nd))
    [mm,aaa]=min(abs(datenum(dates,'yyyymmdd')-times(nd)));
    data=dates{aaa}
    
    folder=['../meteotsunamis/meteotsunami_',data,'/'];
    
load([folder,'prop_vel_corr_allfreq_',data,'.mat'],'estacions','vvv','phipp','tt')

%%

%--------------------------------------------------------------------------
n=2;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi)+cc,wid,hi-cc];
loc=[lm, bm+(sp+hi2)*(nn-1)+(n-1)*hi+cc,wid,hi-cc];
ax0=subplot('Position',loc);
aux=sl.time>=t_in(nd) & sl.time<=t_out(nd);
plot(sl.time(aux),sl.SL(aux)/100-nanmean(sl.SL(aux)/100))
xlim([t_in(nd) t_out(nd)])
title(ax0,['(',lletres{nn},') ',datestr(t_in(nd),'mmm-dd')])
ax0.TitleHorizontalAlignment='left';
yl=max(abs(sl.SL(aux)/100));
yl=0.7;
ylabel('SL [m]')
ylim([-yl-0.05 yl+0.05])
% datetick
ax=gca;
xl=xlim;
ticks=xl(1):2/24:xl(end);
ax.XTick=ticks;
ax.XTickLabel='';
grid on


%--------------------------------------------------------------------------

pp=[]; pv=[];

% vvv=vvv(1:end-1);
% phis=phipp(1:end-1);

np=1;
vs=vvv;
phis=phipp;




n=1;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi+cc];
loc=[lm, bm+(sp+hi2)*(nn-1)+(n-1)*hi,wid,hi+cc];
ax1=subplot('Position',loc);
hold on


pv=[pv plot(ax1,tt,nanmedian(vs,2),'LineWidth',2)];


aux=not(isnan(nanmedian(vs,2)));
patch(ax1,[tt(aux) flip(tt(aux))],[prctile(vs(aux,:),10,2);flip(prctile(vs(aux,:),90,2))],'b','FaceAlpha',0.1,'LineStyle','none')

ylim([0 50])
ylabel('Prop. speed [m/s]')
grid on
xlim([tt(1) tt(end)])
ax1.YTick=[0 10 20 30 40];
% title(ax1,datestr(times(nd),'dd-mmm-yyyy'))

%xlim([datenum(2021,5,23,7,0,0) datenum(2021,5,24,3,0,0)])
% xlim([times(nd)-wint/2 times(nd)+wint/2])
xlim([t_in(nd) t_out(nd)])

xl=xlim;
ticks=xl(1):2/24:xl(end);
ax1.XTick=ticks;
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
ylim([130 330])

% xlim([times(nd)-wint/2 times(nd)+wint/2])
xlim([t_in(nd) t_out(nd)])
xl=xlim;
ticks=xl(1):2/24:xl(end);
ax2.XTick=ticks;
ax2.XTickLabel=datestr(ticks','HH');
ax2.YTick=[130 170 210 250 290 330];
% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end



% 
% plot(ax1,ww.time,ww.wspeed,'k','LineWidth',1,'Marker','s','MarkerFaceColor','g')
% plot(ax2,ww.time,270-ww.wdir*180/pi,'k','LineWidth',1,'Marker','s','MarkerFaceColor','g')

% per=strcat(num2str(round(limP(2:end-1)')),' - ');
% per=strcat(per,num2str(round(limP(1:end-2)')));

% per=cellstr(strcat(per,' min'));


% legend(ax1,pv,per,'Location','Northeast')
if nn==1
legend(ax2,[pv pp],{'Propagation speed','Propagation direction'},'Location','southwest')
end

nn=nn+1;

end
print(gcf,'-djpeg','-r500',['Fig7.jpg'])
return
%%
N=1;

lm=0.07;
bm=0.04;
sp=0.0;
wid=0.86;
tp=0.03;
sp=0.08

hi2=(1-bm-tp-sp*(N-1))/N;

hi=hi2/2;
cc=0.15/N;

figure('units','centimeters','Position',[2 2 17 7])
nn=1;
lletres=flip({'a','b','c'});

for nd=[6]
    datestr(times(nd))
    [mm,aaa]=min(abs(datenum(dates,'yyyymmdd')-times(nd)));
    data=dates{aaa};
    

load(['prop_vel_corr_allfreq_',data,'.mat'],'estacions','vvv','phipp','tt')

%%
%--------------------------------------------------------------------------
n=2;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi)+cc,wid,hi-cc];
loc=[lm, bm+(sp+hi2)*(nn-1)+(n-1)*hi+cc,wid,hi-cc];
ax0=subplot('Position',loc);
aux=sl.time>=t_in(nd) & sl.time<=t_out(nd);
plot(sl.time(aux),sl.SL(aux)-nanmean(sl.SL(aux)))
xlim([t_in(nd) t_out(nd)])
% title(ax0,['(',lletres{nn},') ',datestr(t_in(nd),'mmm-dd')])
ax0.TitleHorizontalAlignment='left';
yl=max(abs(sl.SL(aux)));
ylabel('SL [m]')
ylim([-yl-0.05 yl+0.05])
% datetick
ax=gca;
xl=xlim;
ticks=xl(1):2/24:xl(end);
ax.XTick=ticks;
ax.XTickLabel='';
grid on


%--------------------------------------------------------------------------

pp=[]; pv=[];

% vvv=vvv(1:end-1);
% phis=phipp(1:end-1);

np=1;
vs=vvv;
phis=phipp;




n=1;
% loc=[lm,(bm+sp*(n-1)+(n-1)*hi),wid,hi+cc];
loc=[lm, bm+(sp+hi2)*(nn-1)+(n-1)*hi,wid,hi+cc];
ax1=subplot('Position',loc);
hold on


pv=[pv plot(ax1,tt,nanmedian(vs,2),'LineWidth',2)];


aux=not(isnan(nanmedian(vs,2)));
patch(ax1,[tt(aux) flip(tt(aux))],[prctile(vs(aux,:),10,2);flip(prctile(vs(aux,:),90,2))],'b','FaceAlpha',0.1,'LineStyle','none')

ylim([0 50])
ylabel('Prop. speed [m/s]')
grid on
xlim([tt(1) tt(end)])
ax1.YTick=[0 10 20 30 40];
% title(ax1,datestr(times(nd),'dd-mmm-yyyy'))

%xlim([datenum(2021,5,23,7,0,0) datenum(2021,5,24,3,0,0)])
% xlim([times(nd)-wint/2 times(nd)+wint/2])
xlim([t_in(nd) t_out(nd)])

xl=xlim;
ticks=xl(1):2/24:xl(end);
ax1.XTick=ticks;
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
ylim([130 330])

% xlim([times(nd)-wint/2 times(nd)+wint/2])
xlim([t_in(nd) t_out(nd)])
xl=xlim;
ticks=xl(1):2/24:xl(end);
ax2.XTick=ticks;
ax2.XTickLabel=datestr(ticks','HH');
ax2.YTick=[130 170 210 250 290 330];
% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end



% 
% plot(ax1,ww.time,ww.wspeed,'k','LineWidth',1,'Marker','s','MarkerFaceColor','g')
% plot(ax2,ww.time,270-ww.wdir*180/pi,'k','LineWidth',1,'Marker','s','MarkerFaceColor','g')

% per=strcat(num2str(round(limP(2:end-1)')),' - ');
% per=strcat(per,num2str(round(limP(1:end-2)')));

% per=cellstr(strcat(per,' min'));


% legend(ax1,pv,per,'Location','Northeast')
if nn==1
legend(ax2,[pv pp],{'Propagation speed','Propagation direction'},'Location','southwest')
end

nn=nn+1;

end
