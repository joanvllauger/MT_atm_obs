% Wavelet spectrum
% Joan Villalonga 30/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     waveletplot_v1(PW,period,time,coi,sig) 
%  
% INPUTS:
%       PW: power wavelete
%       period: periods vector
%       time: time vector
%       coi: valid data cone
%       sig: significance level
%  
%     opcionals
%  
%         Position: [xpos ypos xwidth ywidth] 
%  
%         lag1: fisrt autocorrelation lag. If lag1=0 it si computed here.
%         Default lag1=0.
%  
%         DJ: the spacing between discrete scales. Default is 0.25.
%           A smaller # will give better scale resolution, but be slower to plot.
%  
%         MotherWl: he mother wavelet function. The choices are 'MORLET', 'PAUL', or 'DOG'
%  
%         PARAM = the mother wavelet parameter.
%               For 'MORLET' this is k0 (wavenumber), default is 6.
%               For 'PAUL' this is m (order), default is 4.
%               For 'DOG' this is m (m-th derivative), default is 2.
%  
%         PAD = if set to 1 (default is 1), pad time series with enough zeroes to get
%            N up to the next higher power of 2. This prevents wraparound
%            from the end of the time series to the beginning, and also
%            speeds up the FFT's used to do the wavelet transform.
%            This will not eliminate all edge effects (see COI below).
%  
%         S0 = the smallest scale of the wavelet.  Default is 2*DT.
%   
%         J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%           to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%   
%  
%   OUTPUTS:
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=waveletplot_v2(PW,period,time,coi,sig,varargin)
%--- Default parametes
Position=[10,10,700,300]; yl=[2 120]; color_lim=[-10 3]; f_fig=1; f_colorbar=1;
FSize=16; f_datetime=1; f_sig=1; f_label=1; CM=0;
f_ylabel=1;
yticks=[5 10 20 30 50 100];
y_tick=1;
LW=1.5;
f_log=1;
f_coi=1;
%--- Assiganam les variables opcionals 
if (mod(length(varargin),2)==1), error('Number of imputs has to be even'); 
else
%--- agafam els imparells del varargin que són els que duran els noms
in=1:2:length(varargin); noms=varargin(in);
%--- agafam els parells que son els que tenen els valors dels par?metres
in=2:2:length(varargin); valors=varargin(in); clear in;
end

%--- Comprovam que tots els par?metres siguin correctes.
PosParam={'Position','yl','color_lim','f_fig','f_colorbar','FSize','f_datetime','f_sig','f_label','CM','f_ylabel','LW','f_log','f_coi','y_tick' ,'yticks'};
%--- Comprovam que cada valor dels par?metres s?n del tipo que toca
PosParamType={'double','double','double','double','double','double','double','double','double','char','double','double','double','double','double','double'};
for i=1:length(noms)
    if contains(noms{i},PosParam)==0, error(['The parameter ',noms{i},' is not valid.']); 
    else
        k=find(strcmp(noms{i},PosParam));
        if strcmp(class(valors{i}),PosParamType{k})==0 
        error(['The value of ',noms{i},' must be ',PosParamType{k},' class, not ',class(valors{i}),' class']); 
        else
        eval([noms{i},'=valors{i};']);
        end
    end
end

%--- plot
if logical(f_fig)
figure('Position',Position)
end

pcolor(time,period,log10(PW)); hold on;  shading interp; axis ij;
if CM==0
colormap(gca,'jet');
else
colormap(gca,brewermap([],CM));

end
if f_sig
    hold on
contour(time,period,sig,[-99,1],'k','LineWidth',LW)
end
ax=gca;
if f_coi
patch([time(1),time',time(end)],[max(period),coi,max(period)],'w','EdgeColor','w'); 
end
if f_log
ax.YScale='log';  
end
if y_tick
ax.Layer='top'; 
ax.YLim=yl; 
ax.YTick=yticks;
ax.YGrid='on'; 
ax.YMinorGrid='off';
end
% title('Time Series with gaps')
if f_datetime
datetick('x')
end
if f_label
xlabel('time [hours]')
end
if f_ylabel
ax.YLabel.String='Period [min]';
end
%datetick('x',12)
if logical(f_colorbar)
colorbar; 
end
%ax.YLim=[2,48];
ax.FontSize=FSize;
xlim([time(1) time(end)])
caxis(color_lim);
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = 0.8 ;
ax.XMinorGrid='off';
ax.XGrid='on';
% ax.MinorGridLineStyle = '--';
% ax.MinorGridColor = 'k';
% ax.MinorGridAlpha = 0.5 ;
end
% title(interp)
% text=["Missing data = "+(R*100)+"%"+newline+"Number of gaps = "+N_gaps];
% annotation('textbox', [0.6, 0.22, 0.1, 0.1], 'String', text,'color','k','FontSize',14,'BackgroundColor','w','FaceAlpha',1);
% fig_name=['F5_wavelet_',interp,'_R',num2str(R*100),'_Ngaps',num2str(N_gaps),'.png'];
% saveas(gcf,fig_name,'png')