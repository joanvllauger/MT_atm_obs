addpath '../aux_functions/'
addpath '../m_map/'

% Correm annex1 un senyal amb un loop canviant PARÀMETRES
%
%%
global noise
%--- Velocitat
vel_vect=[15:5:40];
deg=30;
deg=deg*pi/180;

%--- Renou 
noises=[0:0.1:0.5];

%--- Periodes
periods=[5 10 15 20 30 60];

%--- Interval de mostrig
dt=60;


%--- flag plot
fplot=1;
%%
err_cell=cell(length(vel_vect),length(noises));
errp_cell=cell(length(vel_vect),length(noises));
bias_cell=cell(length(vel_vect),length(noises));
biasp_cell=cell(length(vel_vect),length(noises));
%%
for n=1:length(vel_vect)
    for k=1:length(noises)
    vel=vel_vect(n)
    noise=noises(k)
    sub_SI1_234_una_senyal
    
    err_cell{n,k}=err;
    errp_cell{n,k}=errp;
    bias_cell{n,k}=bias;
    biasp_cell{n,k}=biasp;
    end
end



% save(['Annex1_vel_renou_DT',num2str(dt)],'noises','vel_vect','err_cell','errp_cell','bias_cell','biasp_cell')

%% Dibuixam
% load(['Annex1_vel_renou_DT60'],'noises','vel_vect','err_cell','errp_cell','bias_cell','biasp_cell')

nn=3
figure('Units','centimeters','Position',[2 2 15 10])
lm=0.08;
bm=0.06;
vvs=0.05;
hs=0.00;
tp=0.03;
rm=0.03;
hi=(1-vvs*1-bm-tp)/2;
wi=(1-lm-rm);

ticks=[1:length(vel_vect)];

n=1
k=1;
pos=[lm+(wi+hs)*(k-1) 1-tp-hi-(hi+vvs)*(n-1) wi hi];
ax1=axes('Position',pos);

a=boxplot(cell2mat(err_cell(:,nn))','positions',ticks-0.12,'widths',0.18,'colors','b');
hold on
b=boxplot(cell2mat(bias_cell(:,nn))','positions',ticks+0.12,'widths',0.18,'colors','r');
ax=gca;
ylim([-2 2])
ax1.YTick=[-2:1:2]
grid on
ax1.XTick=ticks;
ax1.XTickLabel='';
% ax1.Color='none';
ylabel('m/s')
% legend('RMSE','BIAS')
h=findall(gca,'Tag','Box');
m=findall(gca,'Tag','Median');
o=findall(gca,'Tag','Outliers');
for j=1:6
  x = patch(get(h(j),'XData'),get(h(j),'YData'), 'r','FaceAlpha',.5);
  set(m(j),'color','k')  
  % MODIFICATION TO SET THE LINE COLOR TO PATCH COLOR:
    h(length(h)-j+1).Color = x.FaceColor; % reordered to match
end
  
for j=7:length(h)
  x = patch(get(h(j),'XData'),get(h(j),'YData'), 'b','FaceAlpha',.5);
  set(m(j),'color','k')
  set(o(j),'color','b','Marker','+','MarkerEdgeColor','b')
    % MODIFICATION TO SET THE LINE COLOR TO PATCH COLOR:
    h(length(h)-j+1).Color = x.FaceColor; % reordered to match
end

hLegend = legend(h([1 7]), {'RSME','BIAS'});
set(hLegend,...
    'Position',[0.0887712878665105 0.545604485338944 0.146384478035848 0.0859788337712565]);
title('(a) RMSE and BIAS of the speed estimates')
% Among the children of the legend, find the line elements
% hChildren = findall(get(hLegend,'Children'), 'Type','Line')
% % Set the horizontal lines to the right colors
% set(hChildren(6),'Color',[1 0 0])
% set(hChildren(4),'Color',[0 0 1])
% set(hChildren(2),'Color',[0 0.5 0])
ax1.TitleHorizontalAlignment='left';

n=2
k=1;
pos=[lm+(wi+hs)*(k-1) 1-tp-hi-(hi+vvs)*(n-1) wi hi];
ax1=axes('Position',pos);

boxplot(cell2mat(errp_cell(:,nn))','positions',ticks-0.12,'widths',0.18,'colors','b')
hold on
boxplot(cell2mat(biasp_cell(:,nn))','positions',ticks+0.12,'widths',0.18,'colors','r')
ax=gca;
ylim([-3 3])
grid on
ax1.XTick=ticks;
ax1.XTickLabel=num2str(vel_vect');
h=findall(gca,'Tag','Box');
m=findall(gca,'Tag','Median');
o=findall(gca,'Tag','Outliers');
for j=1:6
  x = patch(get(h(j),'XData'),get(h(j),'YData'), 'r','FaceAlpha',.5);
  set(m(j),'color','k')  
  % MODIFICATION TO SET THE LINE COLOR TO PATCH COLOR:
    h(length(h)-j+1).Color = x.FaceColor; % reordered to match
end
  
for j=7:length(h)
  x = patch(get(h(j),'XData'),get(h(j),'YData'), 'b','FaceAlpha',.5);
  set(m(j),'color','k')
  set(o(j),'color','b','Marker','+','MarkerEdgeColor','b')
    % MODIFICATION TO SET THE LINE COLOR TO PATCH COLOR:
    h(length(h)-j+1).Color = x.FaceColor; % reordered to match
end
xlabel('Propagation speed [m/s]')
ylabel('Degrees')
title('(b) RMSE and BIAS of the direction estimates')
ax1.TitleHorizontalAlignment='left';
print(gcf,'-djpeg','-r500',['SI1_4'])
%%
figure
boxplot(cell2mat(bias_cell(:,2)'))
ax=gca;
ax.XTickLabel=num2str(vel_vect');


figure
boxplot(cell2mat(errp_cell(:,2)'))
ax=gca;
ax.XTickLabel=num2str(vel_vect');

figure
boxplot(cell2mat(biasp_cell(:,2)'))
ax=gca;
ax.XTickLabel=num2str(vel_vect');