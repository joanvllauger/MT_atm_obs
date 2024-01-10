addpath '../aux_functions/'
% 1--- Legim les dades de nivell del mar de Ciutadalla
% load('waveheigth_ciutadella_VENOM.mat','time','SL')
% time_sl=time(31:end);
% sl=SL(31:end);
% 
T=readtable('../List_of_events.xlsx');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));
add={'','','','','','a','b','','',''};

load('../data/ciutadella_SL_AtmPres.mat','time','SL','Patm','lon','lat','mareografs')

%% --- coses del plot
SL=SL/100;

lm=0.09;
bm=0.08;
sp=0.01;
tp=0.06;
wid=0.83;
hi=(1-3*sp-bm-tp)/4;
box=[0.04 0.04];
f_size=12;


%% Dibuix
bm=0.035;
tp=0.024;
lm=0.06;
rm=0.09;

is=0.005;

hs=0.005;
vs=0.03-is;

f_size=10;

N=5;
wi=(1-lm-rm-hs)/2;
hi=(1-tp-bm-vs*(N-1))/N;


h=figure('Units','centimeters','Position',[2 2 33 N*4-1]); 
for ne=1:N%length(times)
    disp(datestr(times(ne)))
auxs=time>=t_in(ne) & time<t_out(ne);

sl=SL(auxs); timesl=time(auxs);
tini=t_in(ne); tfi=t_out(ne);
y=rellenar_huecos(timesl,sl,timesl(end),'lineal');
nans=isnan(sl);
sl=F2_filt_simple(timesl,y,2,2*60);
sl(nans)=NaN;
% DT=floor((tfi-tini)*24/4)/24;
DT=5/24;


ticks=[round((tini)*24+2)/24:DT:tfi-DT/4];
% Dibuixam SL
n=3;
axes('position',[lm 1-tp-(hi+vs)*(ne-1)-hi wi hi/2-is])

plot(timesl,sl,'LineWidth',1);%,'o'
hold on
ax=gca;
ax.FontSize=f_size;
% xlabel('Time UTC')
% ylabel('SL [m]')
datetick('x')
xlim([tini,tfi])
%ylim([-0.7 0.7])
grid on
%title(mareograf,'Interpreter','none')

%--- Ylim 
mx=max(abs(sl),[],'all');
if mx<0.7
ylim([-0.7 0.7])
ytick=[-0.5 0 0.5];
else
    ylim([-mx-0.05 mx+0.05])
end
ax.YTick=ytick;

%---
ax.XTick=ticks;
ax.XTickLabel=datestr(ticks,'hh');
ax.YMinorTick='off';
ax.YMinorGrid='off';
ax.XMinorTick='off';
ax.XMinorGrid='on';
ax.GridAlpha = 0.8 ;
ax.MinorGridColor = 'k';
ax.MinorGridAlpha = 0.8 ;
ylab=annotation('textbox',[lm*0.05 1-tp-(hi+vs)*(ne-1)-hi/15 lm*0.75 hi/5],'String',[datestr(t_in(ne),'mmm-dd'),add{ne}],'Interpreter','none','LineStyle','none');

% ypos=ylab.Position;
% ylab.Position(2)=ypos(2)+0.55;
ylab.FontWeight='bold';

ylab=ylabel('[m]');
ylab.Position(1)=tini-0.06*(tfi-tini);
%--------------------------------------------------------------------------

n=1;
axes('position',[lm+hs+wi 1-tp-(hi+vs)*(ne-1)-hi wi hi/2-is])
% Dibuixam wavelet
[PW,period,coi,sig]=wavelet_v1(y,1);
PW(:,nans)=NaN;
waveletplot_v2(PW,period,timesl,coi,sig,'color_lim',[-4 0.5],'f_sig',0,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0,'LW',0.8)
ax=gca;
ylabel('')
ax.YTick=[5 10 20 50 120];
ax.YAxis.FontSize=f_size-2;
ax.XTick=ticks;
% ax.XTickLabel=[];
ax.XTickLabel=datestr(ticks,'hh');
ax.YAxisLocation='right';
colormap(ax,pmkmp(100,'bonito'));
if ne==1
cb1=colorbar('location','east');
end
%waveletplot_v2(WL{ne}(:,auxt),period,time,WL_coi{ne}(:,auxt),WL_sig{ne}(:,auxt),'color_lim',[-7 0],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0)

% title(mareograf,'Interpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Reetallem el temps en l'interval dessitjat
auxp=(time>=tini & time<tfi);
timep=time(auxp)';
y=rellenar_huecos(timep,Patm(auxp),timep(end),'lineal');
P_filt=F2_filt_simple(time,y,2,2*60);
%---------------------------------------------------------------------------
%

n=4;
axes('position',[lm 1-tp-(hi+vs)*(ne-1)-hi/2 wi hi/2-is])
plot(timep,P_filt,'r','LineWidth',1);%,'o'
hold on
ax=gca;
ax.FontSize=f_size;

%--- Ylim
mx=max(abs(P_filt),[],'all');
if mx<1
ylim([-1.1 1.1])
ytick=[-1:1:1];
else
    ylim([-mx-0.1 mx+0.1])
    ytick=[-2:1:2];
end
ax.YTick=ytick;

% xlabel('Time UTC')
% ylabel('Patm [hPa]')
datetick('x')
xlim([tini,tfi])
grid on
% ylim([-0.7 0.7])
ax=gca;

ylab=ylabel('[hPa]');
ylab.Position(1)=tini-0.06*(tfi-tini);
ax.YMinorTick='off';
ax.YMinorGrid='off';
ax.XMinorTick='off';
ax.XMinorGrid='on';
ax.GridAlpha = 0.8 ;
ax.MinorGridColor = 'k';
ax.MinorGridAlpha = 0.8 ;
ax.XTick=ticks;
ax.XTickLabel=[];

%--------------------------------------------------------------------------
%

n=2;
axes('position',[lm+hs+wi 1-tp-(hi+vs)*(ne-1)-hi/2 wi hi/2-is])
% Dibuixam wavelet
[PW,period,coi,sig]=wavelet_v1(y,1);
% figure
% subplot(211)
% waveletplot_v2(PW,periodd,time,coi,sig,'color_lim',[-5 4],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0)
% subplot(212)
waveletplot_v2(PW,period,timep',coi,sig,'color_lim',[-2.5 2],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0,'LW',0.8)
ax=gca;
if ne==1
cb2=colorbar('location','east');
end
ax.YTick=[5 10 20 50 120];

ax.YAxisLocation='right';
ax.XTick=ticks;
ax.XTickLabel=[];
ax.YAxis.FontSize=f_size-2;
colormap(ax,pmkmp(100,'precip'));
ylabel('')
% title(mareograf,'Interpreter','none')


end
llm=0.03;
hii=0.4;
cb2.Position=[lm+wi*2+hs+llm (1-2*hii-0.1)/2 0.01 hii];
cb1.Position=[lm+wi*2+hs+llm (hii+0.1)+(1-2*hii-0.1)/2 0.01 hii];
cb1.Label.String='log_{10}([m]^2)';
cb2.Label.String='log_{10}([hPa]^2)';
% cb2.Label.Position(2)=cb2.Label.Position(2)+0.8;
% cb1.Label.Position(2)=cb1.Label.Position(2)+0.8;
print(gcf,'-djpeg','-r500','Fig2A')

%% Dibuix
bm=0.035;
tp=0.024;
lm=0.06;
rm=0.09;

is=0.005;

hs=0.005;
vs=0.03-is;

f_size=10;

N=5;
wi=(1-lm-rm-hs)/2;
hi=(1-tp-bm-vs*(N-1))/N;


h=figure('Units','centimeters','Position',[2 2 33 N*4-1]); 
for ne=6:length(times)
    nnne=ne-5;
    disp(datestr(times(ne)))
auxs=time>=t_in(ne) & time<t_out(ne);

sl=SL(auxs); timesl=time(auxs);
tini=t_in(ne); tfi=t_out(ne);
y=rellenar_huecos(timesl,sl,timesl(end),'lineal');
nans=isnan(sl);
sl=F2_filt_simple(timesl,y,2,2*60);
sl(nans)=NaN;
% DT=floor((tfi-tini)*24/4)/24;
DT=5/24;


ticks=[round((tini)*24+2)/24:DT:tfi-DT/4];
% Dibuixam SL
n=3;
axes('position',[lm 1-tp-(hi+vs)*(nnne-1)-hi wi hi/2-is])

plot(timesl,sl,'LineWidth',1);%,'o'
hold on
ax=gca;
ax.FontSize=f_size;
% xlabel('Time UTC')
% ylabel('SL [m]')
datetick('x')
xlim([tini,tfi])
%ylim([-0.7 0.7])
grid on
%title(mareograf,'Interpreter','none')

%--- Ylim 
mx=max(abs(sl),[],'all');
if mx<0.7
ylim([-0.7 0.7])
ytick=[-0.5 0 0.5];
else
    ylim([-mx-0.05 mx+0.05])
end
ax.YTick=ytick;

%---
ax.XTick=ticks;
ax.XTickLabel=datestr(ticks,'hh');
ax.YMinorTick='off';
ax.YMinorGrid='off';
ax.XMinorTick='off';
ax.XMinorGrid='on';
ax.GridAlpha = 0.8 ;
ax.MinorGridColor = 'k';
ax.MinorGridAlpha = 0.8 ;
ylab=annotation('textbox',[lm*0.05 1-tp-(hi+vs)*(nnne-1)-hi/10 lm*0.75 hi/5],'String',[datestr(t_in(ne),'mmm-dd'),add{ne}],'Interpreter','none','LineStyle','none');

% ypos=ylab.Position;
% ylab.Position(2)=ypos(2)+0.55;
ylab.FontWeight='bold';

ylab=ylabel('[m]');
ylab.Position(1)=tini-0.06*(tfi-tini);

%--------------------------------------------------------------------------

n=1;
axes('position',[lm+hs+wi 1-tp-(hi+vs)*(nnne-1)-hi wi hi/2-is])
% Dibuixam wavelet
[PW,period,coi,sig]=wavelet_v1(y,1);
PW(:,nans)=NaN;
waveletplot_v2(PW,period,timesl,coi,sig,'color_lim',[-4 0.5],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0,'LW',0.8)
ax=gca;
ylabel('')
ax.YTick=[5 10 20 50 120];
ax.YAxis.FontSize=f_size-2;
ax.XTick=ticks;
% ax.XTickLabel=[];
ax.XTickLabel=datestr(ticks,'hh');
ax.YAxisLocation='right';
colormap(ax,pmkmp(100,'bonito'));
if nnne==1
cb1=colorbar('location','east');
end
%waveletplot_v2(WL{ne}(:,auxt),period,time,WL_coi{ne}(:,auxt),WL_sig{ne}(:,auxt),'color_lim',[-7 0],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0)

% title(mareograf,'Interpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Reetallem el temps en l'interval dessitjat
auxp=(time>=tini & time<tfi);
timep=time(auxp)';
y=rellenar_huecos(timep,Patm(auxp),timep(end),'lineal');
P_filt=F2_filt_simple(time,y,2,2*60);
%---------------------------------------------------------------------------
%

n=4;
axes('position',[lm 1-tp-(hi+vs)*(nnne-1)-hi/2 wi hi/2-is])
plot(timep,P_filt,'r','LineWidth',1);%,'o'
hold on
ax=gca;
ax.FontSize=f_size;

%--- Ylim
mx=max(abs(P_filt),[],'all');
if mx<1
ylim([-1.1 1.1])
ytick=[-1:1:1];
else
    ylim([-mx-0.1 mx+0.1])
    ytick=[-2:1:2];
end
ax.YTick=ytick;
ylab=ylabel('[hPa]');
ylab.Position(1)=tini-0.06*(tfi-tini);
% xlabel('Time UTC')
% ylabel('Patm [hPa]')
datetick('x')
xlim([tini,tfi])
grid on
% ylim([-0.7 0.7])
ax=gca;


ax.YMinorTick='off';
ax.YMinorGrid='off';
ax.XMinorTick='off';
ax.XMinorGrid='on';
ax.GridAlpha = 0.8 ;
ax.MinorGridColor = 'k';
ax.MinorGridAlpha = 0.8 ;
ax.XTick=ticks;
ax.XTickLabel=[];

%--------------------------------------------------------------------------
%

n=2;
axes('position',[lm+hs+wi 1-tp-(hi+vs)*(nnne-1)-hi/2 wi hi/2-is])
% Dibuixam wavelet
[PW,period,coi,sig]=wavelet_v1(y,1);
% figure
% subplot(211)
% waveletplot_v2(PW,periodd,time,coi,sig,'color_lim',[-5 4],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0)
% subplot(212)
waveletplot_v2(PW,period,timep',coi,sig,'color_lim',[-2.5 2],'f_sig',1,'f_fig',0,'f_colorbar',0,'FSize',f_size,'f_label',0,'LW',0.8)
ax=gca;
if nnne==1
cb2=colorbar('location','east');
end
ax.YTick=[5 10 20 50 120];

ax.YAxisLocation='right';
ax.XTick=ticks;
ax.XTickLabel=[];
ax.YAxis.FontSize=f_size-2;
colormap(ax,pmkmp(100,'precip'));
ylabel('')
% title(mareograf,'Interpreter','none')


end
llm=0.03;
hii=0.4;
cb2.Position=[lm+wi*2+hs+llm (1-2*hii-0.1)/2 0.01 hii];
cb1.Position=[lm+wi*2+hs+llm (hii+0.1)+(1-2*hii-0.1)/2 0.01 hii];
cb1.Label.String='log_{10}([m]^2)';
cb2.Label.String='log_{10}([hPa]^2)';
% cb2.Label.Position(2)=cb2.Label.Position(2)+0.8;
% cb1.Label.Position(2)=cb1.Label.Position(2)+0.8;
print(gcf,'-djpeg','-r500','Fig2B')


return