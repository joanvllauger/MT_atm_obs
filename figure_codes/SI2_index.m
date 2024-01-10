% loop C7
addpath '../aux_functions/'
%%
%--- importam les dades de nivell del mar de Ciutadella
load('../data/ciutadella_SL_long.mat','time','SL')

wlen=60;
% wlen=windays*24*60;

%--- Carregam les dades de ERA5
era=load('../data/sepic_index_vars.mat');
%%
T=readtable('../List_of_events.xlsx');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));


%--- Comencem el vector de temps amb la primera dada que no son NaN
in=find(not(isnan(SL)),1,'first'); fi=find(not(isnan(SL)),1,'last');
SL=SL(in:fi); time=time(in:fi);
SL=QC_SL_v1(SL);
%--- Calculam l'alçada d'ona
[t,wh,t_raw,wh_raw]=wave_height_v1(time,SL,'finestra',wlen);


%%


bm=0.1;
tp=0.052;
lm=0.1;
rm=0.12;


hs=0.03;
vs=0.068;

N=1;
wi=(1-lm-rm);
hi=(1-tp-bm);

%--- Feim nomes un interval de temps
tini=datenum(2012,1,1); tfi=datenum(2022,0,0);
aux=t>=tini&t<=tfi;
t=t(aux); wh=wh(aux);
c=3;
windays=3;

%--- Calculem l'index de Sepic
[wh_filt,S_ind,S_filt,parameters]=FSepic_index_v3(era.time,era.MU,t,wh,'windays',windays,'smoothmax',0);

%% --- Figure 2: scatter serie temporal completa (new)
aux=not(isnan(wh_filt));
tt=t(aux);
vdum1=(S_ind(aux)-nanmean(S_ind))/nanstd(S_ind);
vdum2=wh_filt(aux);

%%
figure('units','centimeters','Position',[2 2 16 12])

th_wh=prctile(vdum2,99.5);
% th_wh=0.7;
th2_wh=prctile(vdum2,95);
aux=vdum2>th_wh;
th=min(vdum1(aux));

%--- Seleccionam els que són part de events
auxe=zeros(length(tt),1);
ttt=zeros(length(times),1);
for nt=1:length(times)
    ind=find(tt>=t_in(nt)&tt<=t_out(nt));
    [a,ind2]=max(vdum2(ind));
    auxe(ind(ind2))=1;
    ttt(nt)=tt(ind(ind2));

end
auxe=logical(auxe);

axes('position',[lm bm wi hi]);

% %--- % over th_wh that is over th;
% oth=vdum2(vdum1>th);
% ratio_mt(nm)=sum(oth>th_wh)/length(oth);
% 
% thv(nm)=th;

scatter(vdum1,vdum2,50,[1 1 1]*0.7,'.')
xlim([-3 4])
ax=gca;
yl=ax.YLim;
hold on
p1=plot(ones(2,1).*th,ylim,'-.b','LineWidth',1.8);
% p2=plot(xlim,ones(2,1).*th_wh,'--b','LineWidth',1.8);
% p3=plot(xlim,ones(2,1).*th2_wh,'--','LineWidth',1.8,'color',[0 0.6 0]);
cc=corrcoef(vdum1,vdum2);
corrv=cc(2,1)
colors=jet(length(ttt))
for nt=1:length(times)
    ind=find(tt>=t_in(nt)&tt<=t_out(nt));
    
    [a,ind2]=max(vdum2(ind));
    auxe(ind(ind2))=1;
    ttt(nt)=tt(ind(ind2));

    scatter(vdum1(ind(ind2)),vdum2(ind(ind2)),35,colors(nt,:),'o','filled','MarkerEdgeColor',colors(nt,:),'MarkerEdgeColor','k')


end

B=robustfit(vdum1,vdum2);
hold on,hp=plot(xlim,xlim*B(2)+B(1),'k');hold off,set(hp,'linewidth',3);
title(['Correlation = ',num2str(cc(2,1),'%1.2f')])

legend([hp,p1],{'Linear regression','index threshold'},'location','northwest')

caxis([1 length(ttt)+1])
grid on

ylim(yl)


cticks=cellstr(datestr(t_in,'mmm-dd'));
cticks{6}=[cticks{6},'a'];
cticks{7}=[cticks{7},'b'];

colormap(jet(length(ttt)))
cb=colorbar;
cb.Ticks=[1:sum(auxe)];
cb.TickLabels=cticks;
cb.Position=[lm+wi+0.01 bm 0.02 hi];

p995=sum(vdum2(vdum1>th)>th_wh)/sum(vdum1>th);
f995=sum(vdum2(vdum1<th)>th_wh)/sum(vdum1<th);

p950=sum(vdum2(vdum1>th)>th2_wh)/sum(vdum1>th);
f950=sum(vdum2(vdum1<th)>th2_wh)/sum(vdum1<th);

fprintf('Percetile 99.5 -> MT: %2.2f fals negatiu: %2.2f \n',p995*100,f995*100)
fprintf('Percetile 95 -> MT: %2.2f fals negatiu: %2.2f \n',p950*100,f950*100)



xlabel('Normalized Index','fontweight','bold','fontsize',12)


ylabel({'Filtered SLO Amplitude [m]'},'fontsize',12,'fontweight','bold')


print(gcf,'-djpeg','-r500',['SI2'])




return
