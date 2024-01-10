addpath '../aux_functions/'
%%
T=readtable('../List_of_events_long.xlsx');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));

load('../data/ciutadella_spectral_data.mat')

% wavelet parameters
k0=12;
ff= (4*pi)/(k0 + sqrt(2 + k0^2));
scale=period/ff;
dj=1/12;
%%

limP=flip(2*2.^([0:5]*7/5));
% limP=flip(2*2.^([0:0.5:5]*5/5));
periods=movmedian(limP,2);
T_vect=periods(2:end);
dt_vect=[2 6 16 30 60 120 180 240 360];
T=10.5;
dT=5;

% Utilitzem la taula d'events escollit per triar els esdeveniments de
% meteotsuinami
f_event=zeros(length(time),1);
for n=1:length(times)
    aux=time>=t_in(n) & time<=t_out(n);
    f_event=f_event | aux;
end



% [f_event]=seleccio_episodis_v4(time,SL/100);
corr_e=zeros(length(dt_vect),length(T_vect));
corr_ne=zeros(length(dt_vect),length(T_vect));
corr_all=zeros(length(dt_vect),length(T_vect));

pol_e=zeros(length(dt_vect),length(T_vect));
pol_ne=zeros(length(dt_vect),length(T_vect));
pol_all=zeros(length(dt_vect),length(T_vect));

corr_12=zeros(length(dt_vect),1);
auxsl=period>limP(end) & period <limP(1);
e_sl=dj*movmean(sum(wsl./scale',1),100);
slfilt=butter_filter_v1(SL,mode(diff(time))*24*60,360,2,'high');
mv=movvar(slfilt,100);
C=mv\e_sl';


figure('units','centimeters','position',[2 2 18 9])
lm=0.08;
bm=0.13;
hi=0.78;
hs=0.045;
rm=0.02;
wi=(1-lm-1*hs-rm)/2;
n=1;


for nd=1:length(dt_vect)

dt=dt_vect(nd);
t=[dt/2:dt/2:length(time)-dt/2];
% t=[1:length(time)];

em_pa=zeros(length(t),length(T_vect));
em_sl=zeros(length(t),length(T_vect));
f_e=f_event(t)==1;

lletres={'a','b'};
%e_sl=sqrt(sum(wsl(aux,:),1));
for nt=1:length(T_vect)
    T=T_vect(nt)

aux=period>limP(nt+1) & period <limP(nt);
% aux=period>8 & period <14;
e_pa=sum(wpa(aux,:)./scale(aux)',1);
e_pa(nans)=NaN;
e_pa=sqrt(dj*movmean(e_pa,dt)/C);
em_pa(:,nt)=e_pa(t);

e_sl=sum(wsl(aux,:)./scale(aux)',1);
e_sl(nans)=NaN;
e_sl=sqrt(dj*movmean(e_sl,dt)/C);
em_sl(:,nt)=e_sl(t);
nnnans=movmean(nans,dt);
notn=not(nnnans(t));

emm_sl=em_sl(f_e & notn,nt);
emm_pa=em_pa(f_e & notn,nt);
corr=corrcoef(emm_sl(emm_sl>prctile(emm_sl,0.5)),emm_pa(emm_sl>prctile(emm_sl,0.5)));
corr_e(nd,nt)=corr(2,1);
p1=em_pa(f_e & notn,nt)\em_sl(f_e & notn,nt);
pol_e(nd,nt)=p1;

corr=corrcoef(em_sl(not(f_e) & notn,nt),em_pa(not(f_e)& notn,nt));
p2=em_pa(not(f_e)& notn,nt)\em_sl(not(f_e)& notn,nt);
pol_ne(nd,nt)=p2;
corr_ne(nd,nt)=corr(2,1);


corr=corrcoef(em_sl(notn,nt),em_pa(notn,nt));
corr_all(nd,nt)=corr(2,1);

p=polyfit(em_pa(:,nt),em_sl(:,nt),1);

b=em_pa(:,nt)\em_sl(:,nt);


pol_all(nd,nt)=p(1);
if (round(T)==10 & dt==240) | (round(T)==10 & dt==30)
ax=axes('position',[lm+(n-1)*(hs+wi) bm wi hi]);
plot(em_pa(f_e,nt),em_sl(f_e,nt),'s','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',6)
hold on
plot(em_pa(not(f_e),nt),em_sl(not(f_e),nt),'v','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',4)
ylim([0 40])
xlim([0 0.4])
if n==1;
    
end
plt1=plot(xlim,xlim*p1,'-k','LineWidth',2);
plt2=plot(xlim,xlim*p2,'--k','LineWidth',2);
xlabel('Atmospheric pressure STD [hPa]')
ylabel('Sea Level STD [cm]')
%     if p1(2) < 0
%           sign = ' - ';
%       else
%           sign = ' + ';
%     end
p1_leg = [num2str(p1,'% 0.2f'),'x   \rho= ',num2str(corr_e(nd,nt),'%0.2f')]; 
%     if p2(2) < 0
%           sign = ' - ';
%       else
%           sign = ' + ';
%     end
p2_leg = [num2str(p2,'% 0.2f'),'x   \rho= ',num2str(corr_ne(nd,nt),'%0.2f')]; 
leg=legend([plt1 plt2],{p1_leg,p2_leg},'Location','Northwest');
ax=gca;
ax.FontSize=10;
title(sprintf('(%s) Window width = %3.0f minutes',lletres{n},dt))
ax.TitleHorizontalAlignment='left';

if n==2;
 ylabel('')
end
n=n+1;

end

end
corr=corrcoef(em_pa(:,1),em_pa(:,2));
corr_12(nd)=corr(2,1);

end
print(gcf,'-djpeg','-r500',['Fig4'])

%%
FS=10;
yl=[0 180]
markers={'s','o','^','v','d'};
colors=flip({'b','r','y','g','m'})
bm=0.09;
tp=0.052;
lm=0.08;
rm=0.02;


hs=0.02;
vs=0.04;

N=2;
wi=(1-lm-rm-hs*(N-1))/N;
hi=(1-tp-bm-vs*(2-1))/2;

figure('units','centimeters','position',[2 2 18.5 13])
c=1; n=1;
ax1=axes('position',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n) wi hi]);
for nt=1:length(T_vect)
plot(dt_vect/60,corr_e(:,nt),'-','color',colors{nt},'Marker',markers{nt},'MarkerFaceColor',colors{nt},'MarkerEdgeColor','k')
hold on
end
grid on
leg=cellstr([num2str(limP(2:end)','%3.0f'),repmat('-',length(T_vect),1),num2str(limP(1:end-1)','%3.0f'),repmat(' min',length(T_vect),1)]);

% legend(leg,'location','Southeast')
ylim([0.2 1])
ax=gca;
ax.YTick=[0.3:0.1:1];
% xlabel('Averaging window [hours]');
ylabel('Correlation coefficient');
ax.FontSize=FS;
ax.XTickLabels='';
title('Meteotsunami Events')
t=annotation('textbox',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n)+hi 0.05 0.045],'String','(a)');
t.LineStyle='none';
t.FontWeight='bold';


c=2; n=1;
ax2=axes('position',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n) wi hi]);
for nt=1:length(T_vect)
plot(dt_vect/60,corr_ne(:,nt),'-','color',colors{nt},'Marker',markers{nt},'MarkerFaceColor',colors{nt},'MarkerEdgeColor','k')
hold on
end
grid on
% legend(leg,'location','Southeast')
ylim([0.2 1])

% xlim([dt_vect(1)-0.01 dt_vect(end)+50])
ax=gca;
ax.YTick=[0.2:0.1:1];
% xlabel('Averaging window [hours]');
% ylabel('Correlation coefficient');
ax.FontSize=FS;
ax.YTickLabels='';
ax.XTickLabels='';
title('"Calm" periods')
t=annotation('textbox',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n)+hi 0.05 0.045],'String','(b)');
t.LineStyle='none';
t.FontWeight='bold';

c=1; n=2;
ax3=axes('position',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n) wi hi]);for nt=1:length(T_vect)
plot(dt_vect/60,pol_e(:,nt),'-','color',colors{nt},'Marker',markers{nt},'MarkerFaceColor',colors{nt},'MarkerEdgeColor','k')
hold on
end
grid on
% legend(leg,'location','Northwest')
% ylim([0.35 0.85])
ylim(yl)
ax=gca;
ax.YTick=[0:20:180];
xlabel('Averaging window [hours]');
ylabel('Regression coefficient');
ax.FontSize=FS;
t=annotation('textbox',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n)+hi 0.05 0.045],'String','(c)');
t.LineStyle='none';
t.FontWeight='bold';

c=2; n=2;
ax4=axes('position',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n) wi hi]);
for nt=1:length(T_vect)
plot(dt_vect/60,pol_ne(:,nt),'-','color',colors{nt},'Marker',markers{nt},'MarkerFaceColor',colors{nt},'MarkerEdgeColor','k')
hold on
end
grid on
l=legend(leg,'location','Northwest');
l.Position(1)=0.81;
l.Position(2)=0.45;
ylim(yl)
% xlim([dt_vect(1)-0.01 dt_vect(end)+50])
t=annotation('textbox',[lm+(c-1)*(hs+wi) bm+(vs+hi)*(2-n)+hi 0.05 0.045],'String','(d)');
t.LineStyle='none';
t.FontWeight='bold';
ax=gca;
ax.YTick=[0:20:180];
xlabel('Averaging window [hours]');
% ylabel('Regression coefficient');
ax.FontSize=FS;
ax.YTickLabels='';

% ax.XTick=dt_vect;
% ax.XTickLabel=num2str(dt_vect')
print(gcf,'-djpeg','-r500',['SI3'])
return
