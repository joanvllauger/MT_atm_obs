addpath '../wavelet_package/';
addpath '../aux_functions/'
% 1--- Legim les dades de nivell del mar de Ciutadalla
T=readtable('../List_of_events');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));

add={'','','','','','a','b','','',''};


load('../data/ciutadella_SL_AtmPres.mat','time','SL','Patm','lon','lat','mareografs')
time_sl=time;
patm=Patm;
sl=SL/100;
%% Background
if 1
y=rellenar_huecos(time_sl,sl,time_sl(end),'lineal');
[PW,period,coi,sig]=wavelet_v1(y,1); 
BG=mean(PW,2);
p_bg=period;


y=rellenar_huecos(time,patm,time(end),'lineal');
[PW,period,coi,sig]=wavelet_v1(y,1); 
BGp=mean(PW,2);
p_bg=period;
end



%% Dibuix
mat_sl=[];
mat_pa=[];
% DT=5/12;

for ne=1:length(times)
auxs=time_sl>=t_in(ne) & time_sl<t_out(ne);

SL=sl(auxs); timesl=time_sl(auxs);
tini=t_in(ne); tfi=t_out(ne);
% SL=F2_filt_simple(timesl,SL,2,2*60);
DT=floor((tfi-tini)*24/5)/24;
ticks=[tini:DT:tfi];


% Dibuixam wavelet
[PW,period,coi,sig]=wavelet_v1(SL,1);
auxss=(time_sl(auxs)>=times(ne)-DT & time_sl(auxs)<times(ne)+DT);
% wl_sl=mean(PW(:,auxss),2);
sp_sl=mean(PW,2);
p_sl=period;

% [P,f,conf,DOF]=welch_v2(SL,512,1);
% [sp_sl,DOF] = logfilt(0.02,f,P,2);
% p_sl=1./f;


mat_sl=[mat_sl sp_sl];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Reetallem el temps en l'interval dessitjat
auxp=(time>=tini & time<tfi);
timep=time(auxp)';
Patm=patm(auxp);


y=rellenar_huecos(timep,Patm,time(end),'lineal');
% figure
% plot(timep,Patm)
[PW,period,coi,sig]=wavelet_v1(y,1);
auxpp=(time(auxp)>=times(ne)-DT & time(auxp)<times(ne)+DT);
sp_pa=mean(sig(:,auxpp),2);
sp_pa=mean(PW(:,auxpp),2);
p_pa=period;

mat_pa=[mat_pa sp_pa];



% figure
% subplot(121)
% % loglog(pp_sl,wl_sl);
% % hold on
% loglog(p_sl,sp_sl);
% ylim([1e-5 1e0])
% xlim([2 360])
% ax=gca;
% ax.XTick=[5 10 20 40 80 120];
% grid on
% subplot(122)
% loglog(p_pa,sp_pa);
% ylim([1e-3 1e1])
% xlim([2 360])
% ax=gca;
% ax.XTick=[5 10 20 40 80 120];
% grid on
end

%%
bm=0.09;
tp=0.05;
lm=0.07;
rm=0.1;

hs=0.06;
vs=-0.25;

N=11;
wi=(1-lm-rm-hs)/2;
hi=(1-tp-bm-vs*(N-1))/N;

colors=lines(N-1);
fac=7*((hi+vs)/hi)*N;
fac2=7*((-vs)/hi);

h=figure('units','centimeter','position',[2 2 16 13]);
 ax=axes(h,'position',[lm+hs+wi 1-tp-(hi+vs)*(N-1)-hi wi 1-tp-bm]);

for n=1:length(times)+1 
 if n==N

semilogx(p_bg,log10(BG),'k','LineWidth',1.2)
hold on
% semilogx(3*ones(2,1),[16 17],'k')
% semilogx([2.7 3.3],16*ones(2,1),'k')
% semilogx([2.7 3.3],17*ones(2,1),'k')
% % Create textbox
% annotation(h,'textbox',...
%     [0.502655629139073 0.875213849287168 0.0370794701986753 0.0468431771894043],...
%     'String',{'10^1'},...
%     'FitBoxToText','off',...
%     'EdgeColor','none');


ylim([-6 N*2-2])
xlim([2 360])
% ax=gca;
ax.YTick=[-6:2:N*2-3];
ax.XTick=[2 5 10 20 40 70 120 240];
grid on

ylab=strcat('10^{',num2str([-6:2:N*2-3]'));

ax.YTickLabel=strcat(ylab,'}');


ax.YMinorGrid='off';
ax.YGrid='on';
ax.XMinorGrid='off';
box on
ax.Color='none';
% ax.YTickLabel='';

annotation('textbox',[1-rm+0.001 0.22 0.02 0.02],'String','BG','LineStyle','none')
title('(b) Sea level [m^2]')
 else
% axes('position',[lm+hs+wi 1-tp-(hi+vs)*(n-1)-hi wi hi])
semilogx(ax,p_sl,log10(mat_sl(:,n))+(N-n)*2,'color',colors(n,:),'LineWidth',1.2);
hold on
% % ylim([-7 0])
% % xlim([2 360])
% % ax=gca;
% ax.XTick=[2 5 10 20 40 70 120 240];
% 
% grid on
% set(gca,'Visible','off')
% 
% ax.YMinorGrid='off';
% ax.YGrid='on';
% ax.XMinorGrid='off';
% box on
% ax.Color='none';
% ax.YTickLabel=''; ax.XTickLabel='';
annotation('textbox',[1-rm+0.001 0.22+(N-n)*0.07 0.02 0.02],'String',[datestr(t_in(n),'mmm-dd'),add{n}],'LineStyle','none')

 end  
ax.TitleHorizontalAlignment = 'left';

end
xlabel('Period [min]')

% Cercam els pics a nel BG
% loglog(p_bg,BG,'k','LineWidth',1.2)


%
%--------------------------------------------------------------------------


% hs=0.04;
% vs=-0.1;

N=11;
% wi=(1-lm-rm-hs)/2;
hi=(1-tp-bm-vs*(N-1))/N;

colors=lines(N);


fac=8*((hi+vs)/hi)*N;
fac2=8*((-vs)/hi);
ax = axes(h,'position',[lm 1-tp-(hi+vs)*(N-1)-hi wi 1-tp-bm]);

fac=3;

for n=1:length(times)+1
    
 if n==N
%  axes('position',[lm 1-tp-(hi+vs)*(n-1)-hi wi 1-tp-bm])

semilogx(ax,p_bg,log10(BGp)+(N-n)*fac,'color','k','LineWidth',1.2);
% semilogx(3*ones(2,1),[23 25],'k')
% semilogx([2.7 3.3],23*ones(2,1),'k')
% semilogx([2.7 3.3],25*ones(2,1),'k')
% % Create textbox
% annotation(h,'textbox',...
%     [0.05 0.875213849287168 0.0370794701986753 0.0468431771894043],...
%     'String',{'10^2'},...
%     'FitBoxToText','off',...
%     'EdgeColor','none');


ylim([-4 N*fac-1])
xlim([2 360])
% ax=gca;
ax.XTick=[2 5 10 20 40 70 120 240];
ax.YTick=[-3:fac:N*fac-1];
% ax.GridAlpha=0.7;
ylab=strcat('10^{',num2str([-3:fac:N*fac-1]'));

ax.YTickLabel=strcat(ylab,'}');

grid on


ax.YMinorGrid='off';
ax.YGrid='on';
ax.XMinorGrid='off';
% ax.YMinorTick=[-4:1:N*fac-1];
box on
ax.Color='none';

title('(a) Atmospheric pressure [hPa^2]')
% annotation('textbox',[1-rm-0.01 1-tp-(hi+vs)*(n-1)-0.1 0.02 0.02],'String','BG','LineStyle','none')

 else
% axes('position',[lm 1-tp-(hi+vs)*(n-1)-hi wi hi])
semilogx(ax,p_pa,log10(mat_pa(:,n))+(N-n)*fac,'color',colors(n,:),'LineWidth',1.2);
hold on
% ylim([-5 3])
% ylim([1e-6 1e0])
xlim([2 360])
% ax=gca;
ax.XTick=[2 5 10 20 40 70 120 240];
% grid on
% set(gca,'Visible','off')
% 
% ax.YMinorGrid='off';
% ax.YGrid='off';
% ax.XMinorGrid='off';
% box off
% ax.Color='none';
% ax.YTickLabel=''; ax.XTickLabel='';
%annotation('textbox',[1-rm-0.01 1-tp-(hi+vs)*(n-1)-0.1 0.02 0.02],'String',datestr(times(n),'mm/dd'),'LineStyle','none')
end  
ax.TitleHorizontalAlignment = 'left';

end
xlabel('Period [min]')

%%
figure1=gcf;
% Create arrow
annotation(figure1,'arrow',[0.639072847682119 0.639072847682119],...
    [0.169042769857434 0.228289205702648]);

% Create textbox
annotation(figure1,'textbox',...
    [0.769211920529801 0.189409368635436 0.0619139072847681 0.0529531568228096],...
    'String','33',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create arrow
annotation(figure1,'arrow',[0.68046357615894 0.698675496688742],...
    [0.260692464358452 0.217922606924644]);

% Create arrow
annotation(figure1,'arrow',[0.776490066225166 0.735099337748344],...
    [0.219959266802444 0.211812627291242]);

% Create textbox
annotation(figure1,'textbox',...
    [0.648350993377483 0.246435845213848 0.0619139072847681 0.0529531568228095],...
    'String','24',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.510933774834437 0.154786150712829 0.0288013245033117 0.052953156822809],...
    'String',{'5',''},...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create arrow
annotation(figure1,'arrow',[0.539735099337748 0.579470198675497],...
    [0.181262729124236 0.171079429735234]);


% Create textbox
annotation(figure1,'textbox',...
    [0.606960264900662 0.118126272912422 0.0619139072847681 0.0529531568228097],...
    'String',{'10.5'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

print(gcf,'-djpeg','-r500','Fig3')


return
