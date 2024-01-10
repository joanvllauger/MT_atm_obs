%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------

LW=1.2;

addpath '../aux_functions/'
% addpath C:\Users\Usuario\Desktop\VENOM\Codis\Analisi\Article_atm\prop_vel

%-- Prova del càlcul de velocitats
%--------------------------------------------------------------------------
folder=dir('..\meteotsunamis\meteo*');
files={folder(:).name}';
dates=cellfun(@(x) strrep(x,'meteotsunami_',''),files,'UniformOutput',0);
% dates(1)=[];
wint=12/24; % dies

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
% aux=ones(10,1); aux(5)=0;
% aux=logical(aux);
% times=times(aux); t_in=t_in(aux); t_out=t_out(aux);
%%
tp=0.03;
bm=0.055;
lm=0.06;
rm=0.01;
sp=0.035;
hs=0.012;

hi=(1-4*sp-bm-tp)/5;
wi=(1-1*hs-lm-rm)/2;

p1=0.2; p2=(1-p1)/2;

figure('units','centimeters','Position',[2 2 25 30])

lletres={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)'};

nh=0;
nv=5;
MS=5;
for nd=1:length(times)
    datestr(times(nd))
    [mm,aaa]=min(abs(datenum(dates,'yyyymmdd')-floor(times(nd))));
    data=dates{aaa};
    
load(['../meteotsunamis/meteotsunami_',data,'/prop_vel_corr_',data,'.mat'],'estacions','periods','limP','vvv','phipp','tt')


%%

%--------------------------------------------------------------------------
if mod(nd,2)==1
nv=nv-1;
nh=0;
else
    nh=nh+1;
end



loc=[lm+(wi+hs)*nh,bm+(sp+hi)*(nv)+hi*p2*2,wi,hi*p1];
ax0=subplot('Position',loc);
aux=sl.time>=t_in(nd) & sl.time<=t_out(nd);
plot(sl.time(aux),(sl.SL(aux)-nanmean(sl.SL(aux)))/100)
xlim([t_in(nd) t_out(nd)])
if nd==6
    annotation('textbox',[loc(1)-0.01 loc(2)+hi*p1-0.005 0.14 0.025],'string',[lletres{nd},' ',datestr(t_in(nd),'mmm-dd'),'a'],'FontSize',8,'LineStyle','none','FontWeight','bold')

elseif nd==7
    annotation('textbox',[loc(1)-0.01 loc(2)+hi*p1-0.005 0.14 0.025],'string',[lletres{nd},' ',datestr(t_in(nd),'mmm-dd'),'b'],'FontSize',8,'LineStyle','none','FontWeight','bold')

else
annotation('textbox',[loc(1)-0.01 loc(2)+hi*p1-0.005 0.14 0.025],'string',[lletres{nd},' ',datestr(t_in(nd),'mmm-dd')],'FontSize',8,'LineStyle','none','FontWeight','bold')
end
yl=max(abs(sl.SL(aux)));
ylim([-yl-0.05 yl+0.05])

% datetick
ax=gca;

xl=xlim;
ticks=xl(1):5/24:xl(end);
ax.XTick=ticks;
ax.XTickLabel='';
grid on
ylim([-0.72 0.72])
ax0.XMinorGrid='on';
if nh>0
    ax.YTickLabel='';
else
    ylabel('SL [m]')
    ax.YLabel.FontSize=8;
end
%--------------------------------------------------------------------------
colors=jet(length(periods)-1);
pp=[]; pv=[];

vvv=vvv(1:end-1);
phis=phipp(1:end-1);
periods=periods(1:end-1);

for np=1:length(periods)
vs=vvv{np};
phis=phipp{np};



if np==1
% n=2;
loc=[lm+(wi+hs)*nh,bm+(sp+hi)*(nv)+hi*p2,wi,hi*p2];
ax1=subplot('Position',loc);
era=plot(ax1,ww.time,ww.wspeed,'k','LineWidth',0.5,'Marker','s','MarkerFaceColor','g','MarkerSize',MS);

hold on

else
end

pv=[pv plot(ax1,tt,nanmedian(vs,2),'color',colors(np,:),'LineWidth',LW)];


aux=not(isnan(nanmedian(vs,2)));
patch(ax1,[tt(aux) flip(tt(aux))],[prctile(vs(aux,:),10,2);flip(prctile(vs(aux,:),90,2))],colors(np,:),'FaceAlpha',0.1,'LineStyle','none')
ax1.YTick=[20 40];
ylim([0 50])
if nh==0 & np==length(periods)

ax1.YLabel.String='Speed [m/s]';
ax1.YLabel.FontSize=8;
disp(datestr(ax1.YLabel.Position(1)))
aass=ax1.YLabel.Position(1);

ax1.YLabel.Position(1)=aass-0.025*(t_out(nd)-t_in(nd));
disp(datestr(ax1.YLabel.Position(1)))
end
if nh>0
    ax1.YTickLabel='';
end
grid on
xlim([tt(1) tt(end)])
% title(ax1,datestr(times(nd),'dd-mmm-yyyy'))

%xlim([datenum(2021,5,23,7,0,0) datenum(2021,5,24,3,0,0)])
% xlim([times(nd)-wint/2 times(nd)+wint/2])
xlim([t_in(nd) t_out(nd)])
ax=gca;
xl=xlim;
ticks=xl(1):5/24:xl(end);
ax.XTick=ticks;
ax.XTickLabel='';
cs=cos(phis); sn=sin(phis);
mphi=atan(sum(sn,2,'omitnan')./sum(cs,2,'omitnan'));
ax.XMinorGrid='on';
phis_de=270-phis*180/pi; mphi=270-mphi*180/pi;

aux=phis_de<=0;
phis_de(aux)=360+phis_de(aux);



if np==1
n=1;
loc=[lm+(wi+hs)*nh,bm+(sp+hi)*(nv),wi,hi*p2];
ax2=subplot('Position',loc);
plot(ax2,ww.time,270-ww.wdir*180/pi,'k','LineWidth',0.5,'Marker','s','MarkerFaceColor','g','MarkerSize',MS)

hold on
vs(:,7)=NaN*ones(length(tt),1);

else 
end
% p0=plot(tt,phis_de,'color',[0.5 0.5 0.5]);
% hold on
% p1=plot(tt,prctile(phis_de,25,2),'-.','LineWidth',1,'color','b');
% hold on
% p2=plot(tt,prctile(phis_de,75,2),'--','LineWidth',1,'color','b');
% hold on
pp=[pp plot(tt,nanmedian(phis_de,2),'color',colors(np,:),'LineWidth',LW)];
aux=not(isnan(nanmedian(phis_de,2)));
patch([tt(aux) flip(tt(aux))],[prctile(phis_de(aux,:),10,2);flip(prctile(phis_de(aux,:),90,2))],colors(np,:),'FaceAlpha',0.1,'LineStyle','none')
% p4=plot(tt,ones(length(tt),length(c)).*phi,'--g','LineWidth',1.5);
% for n=1:length(c)
% p4=plot(tt,ones(length(tt),1).*phi(n),'--','LineWidth',1.5,'color',colors(n,:));
% end
% p4=plot(time(auxt),alph_p,'-s','LineWidth',1,'MarkerEdgeColor','k','color','g')

% datetick

if nh==0
ylabel('Dir. [º]')
ax2.YLabel.FontSize=8;
else 
    ax2.YTickLabel='';
end


grid on
ylim([130 300])

% xlim([times(nd)-wint/2 times(nd)+wint/2])
xlim([t_in(nd) t_out(nd)])
ax=gca;
xl=xlim;
ticks=xl(1):5/24:xl(end);
ax.XTick=ticks;
ax.XTickLabel=datestr(ticks','HH');
ax.XMinorGrid='on';
% for nv=1:length(vert_lines)
%     plot(ones(2,1)*vert_lines(nv),ylim,'--k','LineWidth',2)
% end



end


per=strcat(num2str(round(limP(2:end-1)')),' - ');
per=strcat(per,num2str(round(limP(1:end-2)')));

per=cellstr(strcat(per,' min'));


% 
% legend(ax2,pp,per,'Location','southwest')


end
per{end+1}='ERA5 wind at 500 hPa';
leg=legend(ax1,[pv era],per,'Location','Northeast','Orientation','horizontal');
wid=0.35
leg.Position=[(1-wid)/2 bm-0.045 wid 0.015];
print(gcf,'-djpeg','-r500',['SI5.jpg'])

