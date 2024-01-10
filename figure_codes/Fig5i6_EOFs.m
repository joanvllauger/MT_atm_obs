addpath '../m_map/';
addpath '../aux_functions/'
% close all; clear all;
% dibuixam les correlacions
NM=11;

load(['../data/corr_rissagues_1min_allfreq_12h.mat'])

T=readtable('../List_of_events.xlsx');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));
wh_max=datenum(table2array(T(:,4)));

%--- Creem el cevtor de P amb totes les rissagues concatenades
ntt=1;

%Paux=P_t{ntt};
Paux=cat(1,P_t{:});

estacions=mareografs;

add={'','','','','','a','b','','',''};

rep=100;
%% Es fa la matriu de correlació
cr=zeros(length(estacions),length(estacions));
lg=zeros(length(estacions),length(estacions));
for ne=1:length(estacions)
    for nei=1:length(estacions)
        auxc=zeros(length(times),1);
        auxl=zeros(length(times),1);
        for nt=1:length(times)
            if isempty(corr_max_cell{nt})
                auxc(nt)=NaN;
                auxl(nt)=NaN;
            else
            auxc(nt)=corr_max_cell{nt}(ne,nei);
            auxl(nt)=corr_mlag_cell{nt}(ne,nei);
            end
        end
        cr(ne,nei)=mean(auxc,'omitnan');
        lg(ne,nei)=mean(auxl,'omitnan');
    end
end

%% Eliminarem totes les files i columnes que tenguin algun NaN
C=cr;
L=lg;

%% Ordernar
[lon,ind]=sort(lon); lat=lat(ind);
estacions=estacions(ind);

veps_all=zeros(length(estacions),length(estacions));
vexp_all=zeros(length(estacions),1);
C=C(ind,ind);
C(isnan(C))=0;

aux=sum(isnan(C),1)==length(estacions);
C=C(not(aux),not(aux)); 
%--- Calculam els valors propis
[veps_all(not(aux),not(aux)),vaps]=eigs(C,length(estacions));
vexp_all(not(aux))=100*diag(vaps)./sum(diag(C));
veps_all(veps_all==0)=NaN;
vexp_all(vexp_all==0)=NaN;

ne=sum(not(aux));
tlen=12*60*length(times);
p=randn(tlen,ne,rep);%.*randn(1,ne);


%--- Calculam els valors propis
vexp=zeros(ne,rep);
for n=1:rep
C=corrcoef(p(:,:,n));
[veps,vaps]=eigs(C,ne);
vexp(:,n)=100*diag(vaps)./sum(diag(C));

end
 vexp=sort(vexp,2,'ascend');
 disp(strcat('Max V.Exp = ',num2str(vexp(1,95))))
 
 siglev=vexp(:,95);

%--- Calculam els valors propis
% [veps,vaps]=eigs(C,ne);
% vexp=100*diag(vaps)./sum(diag(C));
% 
% disp(['Max V.Exp = ',num2str(vexp(1))])

%% EOF de l'esdeveniment  amb corr màxima

veps_event=zeros(length(estacions),length(estacions),length(times));
vexp_event=zeros(length(estacions),length(times));
%% Ordernar
 
for ntt=1:length(times)

C=corr_max_cell{ntt};
C=C(ind,ind);

%--- Eliminem les finales que no volem
aux=sum(isnan(C),1)==length(estacions);
C=C(not(aux),not(aux)); 

C(isnan(C))=0; 
%--- Calculam els valors propis
[veps,vaps]=eigs(C,sum(not(aux)));
veps_event(not(aux),not(aux),ntt)=veps;

vexp_event(not(aux),ntt)=100*diag(vaps)./sum(diag(C));

[vexp_event(:,ntt),id]=sort(vexp_event(:,ntt),'descend');
veps_event(:,:,ntt)=veps_event(:,id,ntt);

ne=sum(not(aux));
tlen=12*60;
p=randn(tlen,ne,rep);%.*randn(1,ne);

%--- Calculam els valors propis
vexp=zeros(ne,rep);
for n=1:rep
C=corrcoef(p(:,:,n));
[veps,vaps]=eigs(C,ne);
vexp(:,n)=100*diag(vaps)./sum(diag(C));

end
 vexp=sort(vexp,2,'ascend');
 disp(strcat('Max V.Exp = ',num2str(vexp(1,95))))
 
 siglev=zeros(length(aux),1);
 siglev(not(aux))=vexp(:,95);

auxx=(vexp_event(:,ntt)<siglev);
vexp_event(auxx,ntt)=NaN;
veps_event(:,auxx,ntt)=NaN;

end
veps_event(veps_event==0)=NaN;
vexp_event(vexp_event==0)=NaN;

%% Correlacionam els modes per cada event amb els modes de tots els events
maxmode=NM;
corrs=zeros(maxmode,length(times));
for ntt=1:length(times)
    veps=veps_event(:,:,ntt);
    for nmode=1:maxmode
        corrs(nmode,ntt)=corr_v2(veps_all(:,nmode),veps(:,nmode));
    end
end

%% Obrim les alçades màximes
% T=readtable('../Llista_rissagues_in_out_2014_2021_Ciutadella.xlsx');
% ttt=datenum(table2array(T(:,1)));
% aux=ismember(ttt,times);
% wh_max=table2array(T(:,4));
% wh_max=wh_max(aux);

%% 
FS=10

nh=1;
nv=2;
lm=0.08;
rm=0.08;
bm=0.14;
tm=0.06;
hs=0.012;
vs=0.1;

hi=(1-(nv-1)*vs-tm-bm)/nv;
wi=(1-(nh-1)*hs-lm-rm)/nh;

figure('units','centimeters','Position',[2 2 20 13])



colors=jet(maxmode);
plt=[];
leg={};


leg=cellstr(strcat('M',num2str([1:NM]')));
ssss=cat(2,cellstr(datestr(t_in,'mmm-dd')),cellstr(add)');
sssss=cellfun(@(aa,bb) [aa,bb],ssss(:,1),ssss(:,2),'UniformOutput',0);


loc=[lm bm wi hi]
axes('Position',loc);
ax=gca;
bar(abs(corrs(1:NM,:)'))
grid on
ylabel('Correlation')
ax.XTickLabel=sssss;
%title('Correlació entre el Mode de cada rissaga amb el mode mitja')
title('(b) Correlation between event mode and mean mode')
ax.TitleHorizontalAlignment='left';
xl=ax.XLim;
ax.FontSize=FS;
loc=[lm bm+vs+hi wi hi]
axes('Position',loc);
ax=gca;
bar(vexp_event(1:NM,:)',0.6,'stacked')
% xl=xlim;
grid on
ylabel('Var. Exp. [%]')
ax.XTick=[1:length(times)];
ax.XLim=xl;
ax.YLim=[0 100];
ax.YTick=[0:10:100];
% datestr(t_in,'mmm-dd');

ax.XTickLabel=sssss;%datestr(t_in,'mmm-dd');

% title('Variança explicada per cada Mode (calculats per cada esdeveniment)')
title('(a) Explained variance by every mode')
ax.TitleHorizontalAlignment='left';
ax.FontSize=FS;


loc=ax.Position;
ax1=axes('Position',loc);


plot([1:length(times)],wh_max,'s','color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',1.3)
ax1.XLim=xl;
% ax1.Visible = 'off';
ax1.YAxisLocation='right';
tics=[0.4:0.2:1.4]';
ax1.Visible='on';
ax1.XTickLabel='';
ax1.YTick=tics;
ylim([0.4 1.4])
ax1.XTick=[1:length(times)];

ax1.YGrid='off';
ylabel(ax1,'Maximum SLO height [m]');
% ax1.GridColor='r';
% ax1.GridAlpha=0.5;

ax1.YTickLabel=num2str(tics);
ax1.Color='none';
%[0.741557742879069 0.804 0.0889999988675118 0.11190475892453]
[ll]=legend(ax,leg,'Location','NorthWest','Orientation','horizontal');
ll.ItemTokenSize=[15 18];
ll.Position=[lm 0.01 wi bm/2]
ax1.FontSize=FS;
ax.FontSize=FS;
ll.FontSize=FS;
ll.FontSize
print(gcf,'-djpeg','-r500','Fig5')


%% Feim un mapa

% aux=lat<10;
% lat(aux)=39.611004;
MS=15;

nh=3;
nv=4;
lm=0.01;
rm=0.07;
bm=0.01;
tm=0.07;
hs=0.01;
vs=tm;

hi=(1-(nv-1)*vs-tm-bm)/nv;
wi=(1-(nh-1)*hs-lm-rm)/nh;

nt=3;
figure('units','centimeters','Position',[2 2 18 18])
n=1;
k=1;

lletres={'a','b','c','d','e','f','g','h','i','j','k','l'};
signes=ones(1,12); signes(2)=-1; signes(9)=1;
signes(8)=-1;
signes(12)=-1;

for nt=[3 4 9 13]%:length(times)
    
%--- Retallam en lon i lat
latlim=[38.5 40.2];
lonlim=[1 4.5];
%--- Dibuixam
%-- nombre del mode
if nt>length(times)
veps=veps_all;
vaps=vexp_all;
data_str='All events mean EOF'
else
    veps=veps_event(:,:,nt);
    vaps=vexp_event(:,nt);
    data_str=[datestr(t_in(nt),'mmm-dd'),' meteotsunami event']
end
for nmode=1:3
    loc=[lm+(nmode-1)*(wi+hs) bm+(nv-n)*(hi+vs) wi hi];
axes('Position',loc);
m_proj('mercator','long',lonlim,'lat',latlim,'rect');
% m_gshhs_i('patch',[.5 .6 .5]);
m_gshhs_i('color','k','LineWidth',1)
hold on
aux=not(isnan(veps(:,nmode)));
m_scatter(abs(lon(aux)),lat(aux),MS,veps(aux,nmode)*signes(k),'filled','marker','o','MarkerEdgeColor','k')
colormap('jet')
% colorbar
title(['(',lletres{k},') Mode = ',num2str(nmode),': Var. Expl.: ',num2str(vaps(nmode),'%2.2f'),'%'])

m_grid('linestyle','none','linewidth',1,'tickdir','in','yticklabels',[],...
              'xticklabels',[],'xaxisloc','top','fontsize',12);      
caxis([-0.2 0.2])
k=k+1;
end
cb=colorbar('eastoutside');
cb.Position=[(hs+wi)*nh+lm bm+(nv-n)*(hi+vs) (rm)/3 hi];


str=data_str;

t=annotation('textbox',[lm bm+(vs+hi)*(nv-n)+hi+0.04 0.5 0.03],'String',str);
t.FontSize=11;
t.FontWeight='bold';
t.LineStyle='none';
l=annotation('line',[lm lm+0.4],[bm+(vs+hi)*(nv-n)+hi+0.033 bm+(vs+hi)*(nv-n)+hi+0.033]);


n=n+1;
end
print(gcf,'-djpeg','-r500','Fig6')

