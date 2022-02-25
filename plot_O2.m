
function plot_O2(p_point,sonde_index,span_days,Sonde,Model,Counts,Range,Time,Options,Temperature,Format,Alpha,cloud_SDm,HSRL,Data,SNRm,cloud_SDm_above,N_wv,N_wv0,N_wvm,N_wv0m,AbsHumm,AbsHum0m)

%==Fitted lapse rate and surface temperature
figure(728590)
subplot(2,1,1)
plot(Time.thr,permute(Temperature.L_fit_sm_test(1,:,end)*1000,[3 2 1]))
hold on
plot(Time.ts([p_point(1) p_point(end)])/60/60,[-10 -4],'k','linewidth',2)
hold off
yline(-6.5);
grid on
title('Fitted lapse rate')
title(sprintf('Fitted lapse rate\n %s',span_days(1)))
xlabel('Time (UTC hr)')
ylabel('Lapse rate (K/km)')
subplot(2,1,2)
plot(Time.thr,Model.T(1,:) - permute(Temperature.Ts_fit(1,:,:),[3 2 1]))
hold on
plot(Time.ts([p_point(1) p_point(end)])/60/60,[-20 20])
hold off
grid on
title('Surface temp - fitted surface temp')
xlabel('Time (UTC hr)')
ylabel('\deltaT (K)')



%==Background Counts
figure(9985)
%semilogy(datenum(Time.date_ts),Counts.bg_o2off,'-',datenum(Time.date_ts),Counts.bg_o2on,datenum(Time.date_ts),Counts.bg_o2off_mol,datenum(Time.date_ts),Counts.bg_o2on_mol,datenum(Time.date_ts([p_point(1) p_point(end)])),[1 100])
semilogy(datenum(Time.date_ts),Counts.bg_o2off,'-',datenum(Time.date_ts),Counts.bg_o2on,datenum(Time.date_ts),Counts.bg_o2off_mol,datenum(Time.date_ts),Counts.bg_o2on_mol,datenum(Time.date_ts),Counts.bg_wvon,datenum(Time.date_ts),Counts.bg_wvoff,datenum(Time.date_ts([p_point(1) p_point(end)])),[1 100])
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
% xline(738397.94)
% xline(738398.17)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title('Background')
xlabel('Time UTC hr')
ylabel('photons')
legend('Offline','Online','Offline molecular','Online molecular','Online WV','Offline WV')
grid on

%=Diff Counts and sonde difference
% figure(99496)
% cla reset
% line(diff(Counts.o2off(:,p_point(1)))/Range.rangeBin,Range.rkm(1:end-1),'Color','r')
% line(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin./Counts.o2off(1:end-1,p_point(1))*100,Range.rkm(1:end-1),'Color','b')
% %line(diag(dBSRSG(:,p_point)./BSR(:,p_point)),Range.rkm,'Color','b')
% %line(diag(dBSRlo(:,p_point)),Range.rkm(1:end-1))
% legend('dOFF/dr')
% xlim([-2 2])
% xlabel('dOFF/dr')
% ylabel('Range (km)')
% ax1 = gca; % current axes
% ax1.XColor = 'r';
% ax1.YColor = 'r';
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% %line(diag(alpha_0(:,p_point)-Sonde.absorption_sonde{sonde_index}),Range.rkm,'Parent',ax2,'Color','k')
% line(diag(Alpha.alpha_0_filt(:,p_point)-Sonde.absorption_sonde{sonde_index}),Range.rkm,'Parent',ax2,'Color','g')
% %line((diag(alpha_0_filt(1:end-1,p_point)))./(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin./Counts.o2off(1:end-1,p_point(1))*100)-Sonde.absorption_sonde{sonde_index}(1:end-1),Range.rkm(1:end-1),'Parent',ax2,'Color','k')
% %line((diag(alpha_0_filt(1:end-1,p_point)))./(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin)-Sonde.absorption_sonde{sonde_index}(1:end-1),Range.rkm(1:end-1),'Parent',ax2,'Color','k')
% grid on
% xlim([-.5 .5]*10^-3)
% xlabel('\alpha_0 - \alpha_{sonde}')
% ylabel('Range (km)')
% legend('\alpha_0 - \alpha_{sonde}')


%sonde Diff
% % figure(883468)
% % hold on
% % for ii=1:length(Sonde.sonde_ind(1,:))
% %     plot(Sonde.T_sonde(:,ii)-diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,ii))),Range.rkm)
% % end
% % hold off
% % grid on
% % xline(1)
% % xline(-1)
% % title('Sonde - MPD temperature')
% % ylabel('Range')
% % xlabel('\Delta T Sonde-MPD')
% % 
% % figure(883470)
% % hold on
% % for ii=1:length(Sonde.sonde_ind(1,:))
% %     plot(Sonde.T_sonde(:,ii)-diag(Model.T(:,Sonde.sonde_ind(:,ii))),Range.rkm)
% % end
% % hold off
% % grid on
% % xline(1)
% % xline(-1)
% % title('Sonde - lapse')
% % ylabel('Range')
% % xlabel('\Delta T Sonde-lapse')

figure(883469)
%[~,sondeCut]=min(abs(Range.rm-250));
%[~,sondeCut]=min(abs(Range.rm-1000));
tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:Range.i_range
        for j = -60:60
            %tempProb(i,j) = tempProb(i,j)+nnz((diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(240+i-1)+1)&(diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(240+i-1)-0)& (Sonde.T_sonde(sondeCut:end,ii)<(240+j-1)+1)&(Sonde.T_sonde(sondeCut:end,ii)>=(240+j-1)-0));
           tempDiff = Sonde.T_sonde(i,ii)-Temperature.T_finalm(i,Sonde.sonde_ind(i,ii));
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
colormap(flipud(hot))%colormap hot
%set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
%plot(239:302,239:302,'k','LineWidth',2)
%plot(240:303,237:300,'k--','LineWidth',2)
%plot(237:300,240:303,'k--','LineWidth',2)
xline(0)
hold off
%xlim([239 302])
%ylim([239 302])
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\Delta T Sonde-MPD (^oC)')
ylabel('Range')
grid on


%=Water Vapor profile
figure(99499)
plot(diag(Model.WV(:,p_point)),Range.rm)
hold on
plot(Sonde.WV_sonde(:,sonde_index),Range.rm)
hold off
%xlim([0 10^22])
legend('Model','sonde')
xlabel('WV (molec/m^3')
ylabel('Range (m)')
legend('VW')
xlabel('Absolute Humidity g/m^3')

%=Water Vapor profile
figure(99410)
plot(diag(N_wv0m(:,p_point)),Range.rkm,'Linewidth',2)
hold on
plot(diag(N_wvm(:,p_point)),Range.rkm,'--','Linewidth',2)
%plot(diag(N_wv0(:,p_point)),Range.rkm)
plot(Sonde.WV_sonde(:,sonde_index),Range.rkm,'.-')
hold off
xlim([1e23 3e23])
ylim([0 4])
grid on
legend('VW0','VW puturbative','raw','Radiosonde')
legend('VW zero order','VW puturbative','Radiosonde')
xlabel('WV (molec/m^3)')
ylabel('Range (km)')

%=Temperature profile
figure(884)
%%%plot(diag(Model.T(:,p_point)),Range.rkm)
%plot(T(:,p_point)+2,Range.rkm)
%plot(T(:,p_point)-2,Range.rkm)
%plot(Ts_fit(1,p_point,end),0,'+')
temp = diag(Temperature.T_final_test(:,p_point)).*diag(SNRm(:,p_point))-273.13;
temp(temp<=0)=nan;
plot(temp,Range.rkm)
hold on
plot(diag(Temperature.T_finalm(:,p_point))-273.13,Range.rkm,'linewidth',2)
%hold on
%Temperature.exclusion(Temperature.exclusion==0)=NaN;
%plot(diag(T_final_test(:,p_point)).*diag(exclusion(:,p_point,end)),Range.rkm,'*')

%plot(Model.Ts(p_point(1)),Range.rkm(1),'+')
%hold on
%plot(Model.T(:,p_point(1))-273.13,Range.rkm,'-','LineWidth',2)
plot(Sonde.T_sonde(:,sonde_index)-273.13,Range.rm/1000,'.-')
%%plot(Sonde.T_sonde(:,sonde_index)+2,Range.rm/1000,'--k')
%%plot(Sonde.T_sonde(:,sonde_index)-2,Range.rm/1000,'--k')
%plot(diag(L_fit_sm_test(:,p_point,end) .* Range.rm + Ts_fit(1,p_point,end)),Range.rkm,'--')
hold off
grid on
xlabel('Temperature (^oC)')
ylabel('Range (km)')
ylim([0 4])
xlim([-10 40])
title(sprintf(['Temp guess and retrieved from pertabative absorption\n' datestr(Time.date_ts(p_point(1)))]))
legend('Retrieved Temperature','Model temperature','Location','northeast')
legend('Retrieved Temperature','Smoothed Temperature','Radiosonde Temperature','Location','northeast')
%legend('Retrieved Temperature','Radiosonde','Radiosonde+2K','Radiosonde-2K','Location','northeast')
%%%saveas(gcf,[reportfile sprintf('temperature%.0f0%.0f.png',floor(Time.thr(p_point(1))),(Time.thr(p_point(1))-floor(Time.thr(p_point(1))))*60)])


figure(8844)
for ii = 1:30
titer(:,ii) = permute(Temperature.Titer(:,p_point(1),ii),[1 3 2]).*diag(SNRm(:,p_point))-273.13;
end
titer(titer<=0) = nan;
plot(titer,Range.rkm)
grid on
xlabel('Temperature (^oC)')
ylabel('Range (km)')
ylim([0 4])
xlim([-10 40])
title(sprintf(['Temperature iteration\n' datestr(Time.date_ts(p_point(1)))]))



%=Tfit-T
% figure(885)
% plot(diag(Temperature.L_fit_sm_test(:,p_point,end) .* Range.rm + Temperature.Ts_fit(1,p_point,end))-diag(Temperature.T_finalm(:,p_point)),Range.rkm)
% line([0 0],[0 6],'Color','k','LineStyle','--')
% line([-1 -1],[0 6],'Color','k','LineStyle','-.')
% line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
% line([-2 -2],[0 6],'Color','k')
% line([2 2],[0 6],'Color','k')
% xlabel('\Delta T (T_{fit} - T_{retrieved}) (K)')
% ylabel('Range (km)')
% legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','northwest')
% ylim([0 6])
% xlim([-10 10])
% title(sprintf(['Temp fit and DIAL difference\n' datestr(Time.date_ts(p_point(1)))]))
% hold off

%=Tsonde-T
figure(888)
plot(Sonde.T_sonde(:,sonde_index) - diag(Temperature.T_finalm(:,p_point)),Range.rkm)
%plot(diag(Model.T(:,p_point)) - diag(T_finalm(:,p_point)),Range.rkm)
hold on
plot(Sonde.T_sonde(:,sonde_index) - diag(Model.T(:,p_point)),Range.rkm)
plot(diag(Temperature.L_fit_sm_test(:,p_point,end) .* Range.rm + Temperature.Ts_fit(1,p_point,end))-diag(Temperature.T_finalm(:,p_point)),Range.rkm)
hold off
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-1 -1],[0 6],'Color','k','LineStyle','-.')
line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color','k')
line([2 2],[0 6],'Color','k')
xlabel('\Delta T (T_{sonde} - T_{retrieved}) (K)')
ylabel('Range (km)')
ylim([0 5])
xlim([-10 10])
title(sprintf(['Temperature difference\n' datestr(Time.date_ts(p_point(1)))]))
legend('MPD','lapse','fit')
%%%saveas(gcf,[reportfile sprintf('tempDiff%.0f0%.0f.png',floor(Time.thr(p_point(1))),(Time.thr(p_point(1))-floor(Time.thr(p_point(1))))*60)])

figure(888999)
plot((Sonde.WV_sonde(:,sonde_index) - diag(N_wv0(:,p_point)))./Sonde.WV_sonde(:,sonde_index)*100,Range.rkm)
%plot(diag(Model.T(:,p_point)) - diag(T_finalm(:,p_point)),Range.rkm)
hold on
plot((Sonde.WV_sonde(:,sonde_index) - diag(N_wvm(:,p_point)))./Sonde.WV_sonde(:,sonde_index)*100,Range.rkm)
plot((Sonde.WV_sonde(:,sonde_index) - diag(N_wv0m(:,p_point)))./Sonde.WV_sonde(:,sonde_index)*100,Range.rkm)
hold off
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-10 -10],[0 6],'Color','k','LineStyle','-.')
line([10 10],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-20 -20],[0 6],'Color','k')
line([20 20],[0 6],'Color','k')
xlabel('(WV_{sonde} - WV_{retrieved})/WV_{sonde}*100 (%)')
ylabel('Range (km)')
ylim([0 5])
xlim([-100 100])
title(sprintf(['WVdifference\n' datestr(Time.date_ts(p_point(1)))]))
legend('MPD 0','MPD total','MPD0')
%%%saveas(gcf,[reportfile sprintf('tempDiff%.0f0%.0f.png',floor(Time.thr(p_point(1))),(Time.thr(p_point(1))-floor(Time.thr(p_point(1))))*60)])

%=2D temperature
figure(886)
r_max_plot = 4; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
r_min_plot = .5; % Max range to plot [km]
[~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=0; %Set transperency so that cloud mask can be seen
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
colorLimits = [min(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1
colorLimits =[-33 40]+273;
colorLimits =[240 320]-273.13;
%colorLimits = [-30 20];

%colorLimits =[5 30]+273;
colorLimits =[5 30];
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Temperature.T_finalm(1:ind_km_max,:)-273.13,'AlphaData',imAlpha,colorLimits) %Plot
hold on
%sonde lines
% %  plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% % %yline(.6)
for ii = 1:size(Sonde.sonde_ind,2)
    plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
end
% % xline(738397.94)
% % xline(738398.17)
% % xline(738400.65)
% % xline(738401.18)
% xline(738407.17,'LineWidth',2)
% xline(738407.63,'LineWidth',2)
% xline(738408.93,'LineWidth',2)
% xline(738409.72,'LineWidth',2)
%colormap('jet')
colors = colormap; % Set colors to colormap
colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colormap(colors); % Make colormap new colors
set(gca, 'YDir','normal') % Set y axis increasing
set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
xticks(datenum(Time.date_ts(1)):Format.tickSpacing:datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title(sprintf('Temperature retrieval')) %Add title
%xlabel('Time')
ylabel('Range (km)')
title(colorbar,'Temperature (^oC)') % Add title to colorbar
hold off


figure(887)
r_max_plot = 4; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
r_min_plot = .5; % Max range to plot [km]
[~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
% imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
% imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
% imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=0; %Set transperency so that cloud mask can be seen
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen

colorLimits =[0 3e23];
%imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Temperature.T_finalm(1:ind_km_max,:)-273.13,'AlphaData',imAlpha,colorLimits) %Plot
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),N_wv(1:ind_km_max,:)-273.13,colorLimits) %Plot
hold on
%sonde lines
% %  plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% % %yline(.6)
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% % end
colormap('jet')
colors = colormap; % Set colors to colormap
colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colormap(colors); % Make colormap new colors
set(gca, 'YDir','normal') % Set y axis increasing
set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
xticks(datenum(Time.date_ts(1)):Format.tickSpacing:datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title(sprintf('WV number density')) %Add title
xlabel('Time')
ylabel('Range (km)')
title(colorbar,'WV (1/m^3)') % Add title to colorbar
hold off

figure(886645)
r_max_plot = 4; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
r_min_plot = .5; % Max range to plot [km]
[~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=0; %Set transperency so that cloud mask can be seen
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
colorLimits = [min(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1

%colorLimits =[-10 10];
colorLimits =[-2 2]*10^-3;
colorLimits =[-1 1]*10^-3;
%imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Temperature.T_final_testf(:,1:ind_km_max)-Temperature.T_final_testg(:,1:ind_km_max),'AlphaData',imAlpha,colorLimits) %Plot
%%%imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),movmean(abs(Temperature.T_final_testf(1:ind_km_max,:)-Temperature.T_final_testg(1:ind_km_max,:)),5,1),colorLimits) %Plot
%imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),movmean(abs(Alpha.alpha_total_rawf(1:ind_km_max,:)-Alpha.alpha_total_rawg(1:ind_km_max,:)),5,1),colorLimits) %Plot
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),movstd(Alpha.alpha_total_rawf(1:ind_km_max,:)-Alpha.alpha_total_rawg(1:ind_km_max,:),5,0,1),colorLimits) %Plot
hold on
%sonde lines
% %  plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% % %yline(.6)
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% % end
colormap('jet')
colors = colormap; % Set colors to colormap
colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colormap(colors); % Make colormap new colors
colormap(redblue)
set(gca, 'YDir','normal') % Set y axis increasing
%set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
xticks(datenum(Time.date_ts(1)):Format.tickSpacing:datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title(sprintf('Temperature retrieval')) %Add title
%xlabel('Time')
ylabel('Range (km)')
title(colorbar,'Temperature (^oC)') % Add title to colorbar
hold off

%=Absorption profile
figure(4)

gg = cloud_SDm_above.*SNRm;
gg(gg<=0)=nan;
% alpha = diag(Alpha.alpha_0(:,p_point)).*diag(SNRm(:,p_point));
% alpha(alpha==0)=nan;
alpha = Alpha.alpha_0;
alpha(isnan(gg))=nan;
plot(diag(alpha(:,p_point)),Range.rkm)
hold on
plot(diag(Alpha.alpha_totalm(:,p_point)),Range.rkm,'LineWidth',2)
%plot(diag(alpha_total_raw(:,p_point)),Range.rkm)
plot(Sonde.absorption_sonde{sonde_index},Range.rkm,'.-')
%plot(diag(Model.absorption(:,p_point)),Range.rkm)
%plot(diag(alpha_d(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(diag(alpha_0_filt(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(diag(alpha_0_filt_mol(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(alpha_0_filt_decon(:,p_point(1)),Range.rkm,'.-')
%%%plot(alpha_0_filt_corr(:,p_point(1)),Range.rkm,'k.-','LineWidth',2)
%%%plot(Model.absorption(:,p_point(1)),Range.rkm,'-','LineWidth',2)
%plot(alpha_d(:,p_point(1)),Range.rkm(1:end-1),'.-')
%plot(alpha_dModel(:,p_point(1)),Range.rkm(1:end-1),'.-')
hold off
grid on
title(sprintf(['Calculated O_2 absorption coefficient\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')
xlim([-2e-4 6e-4])
ylim([0 4])
%legend('Measured zeroth order absorption','Pertabative absorption','Alpha 0+1+2','Radiosonde model','\alpha model','\alpha_0+\Delta','Smooth \alpha_0')
legend('Measured zeroth order absorption','Smoothed Purtabative absorption','Radiosonde model')
legend('Measured zeroth order absorption','Smoothed Purtabative absorption','Radiosonde Absorption model')
%%%saveas(gcf,[reportfile sprintf('absorption%.0f0%.0f.png',floor(Time.thr(p_point(1))),(Time.thr(p_point(1))-floor(Time.thr(p_point(1))))*60)])

%=Absorption profile
figure(5)
plot(diag(Alpha.alpha_0wv(:,p_point)),Range.rkm);
%plot(Alpha.alpha_0wv(:,p_point(1)),Range.rkm);
hold on
plot(diag(Alpha.alpha_total_rawwv(:,p_point)),Range.rkm);
hold on
hold off
grid on
title(sprintf(['Calculated WV absorption coefficient\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')
xlim([-2e-4 6e-4])
ylim([0 5])
legend('Measured zeroth order absorption','Purtabative absorption')

figure(55)
plot(diag(N_wv0(:,p_point)),Range.rkm);
hold on
plot(N_wv(:,p_point(1)),Range.rkm);
plot(Model.WV(:,p_point(1)),Range.rkm);
plot(Sonde.WV_sonde(:,sonde_index),Range.rkm);
hold off
grid on
title(sprintf(['Calculated WV absorption coefficient\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')
xlim([0 3e23])
ylim([0 5])
legend('Measured zeroth order absorption','Smoothed Purtabative absorption','Radiosonde Absorption model')

% figure(44)
% plot((Model.absorption(:,p_point(1))-Alpha.alpha_totalm(:,p_point(1)))./Model.absorption(:,p_point(1))*100,Range.rkm,'-')
% hold on
% plot((Model.absorption(:,p_point(1))-Alpha.alpha_0(:,p_point(1)))./Model.absorption(:,p_point(1))*100,Range.rkm,'-')
% %%plot((Model.absorption(1:end-1,p_point(1))-alpha_d(:,p_point(1)))./Model.absorption(1:end-1,p_point(1))*100,Range.rkm(1:end-1),'-')
% hold off
% legend('model total diff','model 0 diff','d')
% xlabel('percent diff from model')
% grid on
% xlim([-200 200])


%=BSR
figure(444)
bsr = diag(HSRL.BSR(:,p_point)).*diag(SNRm(:,p_point));
bsr(bsr==0)=nan;
plot(bsr,Range.rkm,'linewidth',2)
hold off
grid on
title(sprintf(['BSR\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlim([0 10])
ylim([0 4])


%=Counts Profile
figure(7)
subplot(2,2,1)
plot(diag(Counts.o2on(:,p_point)),Range.rm)
hold on
plot(diag(Counts.o2off(:,p_point)),Range.rm)
plot(diag(Counts.o2on_mol(:,p_point)),Range.rm)
plot(diag(Counts.o2off_mol(:,p_point)),Range.rm)
plot(Counts.o2on_bgsub(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2off_bgsub(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2on_bgsub_mol(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2off_bgsub_mol(:,p_point(1)),Range.rm_raw_o2,'--')
xlim([0 500])
%view([90 -90])
title('Averaged photon returns')
legend('Online','Offline','Online molecular','Offline molecular','raw online','raw offline','raw online molecular','raw offline molecular')
hold off
grid on
ylabel('Range [m]')
xlabel('Photons')
subplot(2,2,2)
plot(log(Counts.o2on_bgsub(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
hold on
plot(log(Counts.o2off_bgsub(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
plot(log(Counts.o2on_bgsub_mol(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
plot(log(Counts.o2off_bgsub_mol(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('log(Photons*range^2)')
ylabel('Range [m]')
grid on
hold off
subplot(2,2,3)
semilogx((Counts.o2on_bgsub(:,p_point(1))),Range.rm_raw_o2/1000);
hold on
semilogx((Counts.o2off_bgsub(:,p_point(1))),Range.rm_raw_o2/1000);
semilogx((Counts.o2on_bgsub_mol(:,p_point(1))),Range.rm_raw_o2/1000);
semilogx(Counts.o2off_bgsub_mol(:,p_point(1)),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('Photons')
ylabel('Range [m]')
grid on
hold off
subplot(2,2,4)
semilogx(Counts.o2on_bgsub(:,p_point(1))+Counts.bg_o2on(:,p_point(1)),Range.rm_raw_o2/1000);
hold on
semilogx(Counts.o2off_bgsub(:,p_point(1))+Counts.bg_o2off(:,p_point(1)),Range.rm_raw_o2/1000);
semilogx(Counts.o2on_bgsub_mol(:,p_point(1))+Counts.bg_o2on_mol(:,p_point(1)),Range.rm_raw_o2/1000);
semilogx(Counts.o2off_bgsub_mol(:,p_point(1))+Counts.bg_o2off_mol(:,p_point(1)),Range.rm_raw_o2/1000);

% semilogx(Counts.o2on_raw(:,p_point(1)*10),Range.rm_raw_o2/1000,'--');
% semilogx(Counts.o2off_raw(:,p_point(1)*10),Range.rm_raw_o2/1000,'--');
legend('on','off','on mol','off mol')
xlabel('Photons')
ylabel('Range [m]')
grid on
hold off

%=Housekeeping temperature
figure(738732)
subplot(4,1,1)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.InsideCell.Temperature)
legend('Inside Cell')
title('Thermocouples')
ylabel('^o C')
subplot(4,1,2)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.OutsideCell.Temperature  )
ylabel('^o C')
legend('Outside Cell')
subplot(4,1,3)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.TSOA.Temperature  )
ylabel('^o C')
legend('TSOA')
subplot(4,1,4)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.RoomTemp.Temperature  )
ylabel('^o C')
legend('Room Temp')
xlabel('Time UTC hr')

%=Housekeeping temperature
% figure(738732)
% subplot(8,1,1)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.InsideCell.Temperature)
% legend('HSRLOven')
% title('Thermocouples')
% ylabel('^o C')
% subplot(8,1,2)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.HVACReturn.Temperature  )
% ylabel('^o C')
% legend('HVAC return')
% subplot(8,1,3)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.HVACSource.Temperature  )
% ylabel('^o C')
% legend('HVACSource')
% subplot(8,1,4)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.Window.Temperature  )
% ylabel('^o C')
% legend('window')
% subplot(8,1,5)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.WindowHeaterFan.Temperature  )
% ylabel('^o C')
% legend('windowHeaterfan')
% subplot(8,1,6)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.VWEtalonHeatSink.Temperature  )
% ylabel('^o C')
% legend('WVetalon heat sink')
% subplot(8,1,7)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.HSRLEtalonHeatSink.Temperature  )
% ylabel('^o C')
% legend('HSRL etalon heat sink')
% subplot(8,1,8)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.RoomTemp    .Temperature  )
% ylabel('^o C')
% legend('optical table')
% xlabel('Time UTC hr')

%=Laser locking wavelenths
figure(738734)
subplot(2,1,1)
plot(Data.Laser.O2Online.TimeStamp,Data.Laser.O2Online.WaveDiff)
hold on
plot(Data.Laser.O2Online.TimeStamp,Data.Laser.O2Offline.WaveDiff)
hold off
legend('Online','Offline')
title('wavelength Difference')
ylim([-5*10^-4 5*10^-4])
ylabel('(nm)')
xlabel('Time hr')
subplot(2,1,2)
plot(Data.Laser.O2Online.TimeStamp  ,Data.Laser.O2Online.WaveDiff-Data.Laser.O2Offline.WaveDiff  )
title('Online wave diff-offline wave diff')
ylim([-5*10^-4 5*10^-4])
grid on
ylabel('(nm)')
xlabel('Time hr')

%=Laser temperatures
figure(738735)
subplot(4,1,1)
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Online.TemperatureActual)
hold on
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Online.TemperatureDesired)
legend('Actual','Desired')
hold off
grid on
title('Online Laser Temperature')
ylabel('(C)')
subplot(4,1,2)
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Offline.TemperatureActual)
hold on
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Offline.TemperatureDesired)
legend('Actual','Desired')
hold off
grid on
title('Offline Laser Temperature')
ylabel('(C)')
subplot(4,1,3)
% plot(Data.Laser.TWSOA.TimeStamp,Data.Laser.TWSOA.TemperatureActual)
% hold on
% grid on
% plot(Data.Laser.TWSOA.TimeStamp,Data.Laser.TWSOA.TemperatureDesired)
% legend('Actual','Desired')
% hold off
% title('TWSOA Laser Temperature')
% ylabel('(C)')
subplot(4,1,4)
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Etalon.O2Etalon.TemperatureActual)
title('Etalon Temperature')
grid on
ylabel('(C)')
xlabel('Time hr')

%=Laser seed power
figure(778298)
plot(Data.Laser.O2Online.TimeStamp,10.^(Data.Laser.O2Online.SeedPower/10))
hold on
plot(Data.Laser.O2Offline.TimeStamp,10.^(Data.Laser.O2Offline.SeedPower/10))
hold off
%ylim([-5.5 -4.8])
title('Laser Seed power at wavemeter')
legend('Online','Offline')

%=Weather Station
figure(83)
subplot(2,1,1)
plot(datenum(Time.date_ts),Model.Ts-273.13)
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
%title('Weather Station Surface Temperature')
%xlabel('Time')
ylabel('^oC')
xlim([datenum(Time.date_ts(1)) datenum(Time.date_ts(end))])
subplot(2,1,2)
plot(datenum(Time.date_ts),Model.Ps)
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title('Weather Station Suface Pressure')
xlabel('Time')
ylabel('Pressure ATM')
xlim([datenum(Time.date_ts(1)) datenum(Time.date_ts(end))])

%=UPS
if ishandle(648)
    figure(648)
    clf(648)
end
figure(648)
line(Data.UPS.all.TimeStamp,Data.UPS.all.BatteryTimeLeft,'Color','r')
line(Data.UPS.all.TimeStamp,Data.UPS.all.HoursOnBattery,'Color','g')
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('Battery remaining (hrs)')
legend('Battery remaining','Hours on battery','Location','Southwest')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(Data.UPS.all.TimeStamp,Data.UPS.all.UPSTemperature,'Parent',ax2,'Color','k')
ylabel('Battery temperature (C)')

%=2d delta
figure(4489)
imagesc(Time.thr,Range.rm/1000,Counts.delta)
hold on
xline(Time.thr(p_point(1)),'K');
hold off
set(gca, 'YDir','normal')
caxis([-1000 1000]/3)
colorbar
colormap(redblue(64))
title('Counts.delta')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

%=2D counts
figure(589)
subplot(2,2,1)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
colorbar
colormap('hot')
title('o2 on unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 500])
colorbar
colormap('hot')
title('o2 on unaveraged ')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,3)
imagesc(Time.thr,Range.rm,Counts.o2off_noise.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
colorbar
title('o2off unaveraged  * r^{ 2}')
subplot(2,2,4)
imagesc(Time.thr,Range.rm,Counts.o2off_noise)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 500])
colorbar
title('o2off unaveraged')

figure(5899)
subplot(2,1,1)
imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([1 1000])
set(gca,'ColorScale','log')
colorbar
title('o2off unaveraged')
subplot(2,1,2)
imagesc(Time.thr,Range.rm,Counts.o2off_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
set(gca,'ColorScale','log')
caxis([1 1000])
colorbar
title('poisson thin averaged')

%=2D molecular counts
figure(590)
subplot(2,2,1)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise_mol.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
colorbar
colormap('hot')
title('o2 on mol unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 500])
colorbar
colormap('hot')
title('o2 on mol unaveraged')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,3)
imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 5e7])
colorbar
colormap('hot')
title('o2 off mol unaveraged * r^{ 2}')
subplot(2,2,4)
imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 75])
colorbar
colormap('hot')
title('o2 off mol unaveraged')


figure(590)
subplot(2,1,1)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
set(gca,'ColorScale','log')
caxis([0 1e3])
colorbar
colormap('hot')
title('o2 on mol unaveraged')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

subplot(2,1,2)
imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
set(gca,'ColorScale','log')
caxis([0 100])
colorbar
colormap('hot')
title('o2 off mol unaveraged')


figure(591)
subplot(2,1,1)
imagesc(Time.thr,Range.rm/1000,Counts.wvon_noise)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
set(gca,'ColorScale','log')
caxis([0 1e3])
colorbar
colormap('hot')
title('wv on unaveraged')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,1,2)
imagesc(Time.thr,Range.rm,Counts.wvoff_noise)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
set(gca,'ColorScale','log')
caxis([0 1000])
colorbar
colormap('hot')
title('wv off unaveraged')

%=2D absorption
figure(8886)
subplot(3,1,1)
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot));
imAlpha=ones(size(Alpha.alpha_0m(1:ind_km_max,:)));
imAlpha(isnan(Alpha.alpha_0m(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(Alpha.alpha_0m(1:ind_km_max,:),[],'all'), max(Alpha.alpha_0m(1:ind_km_max,:),[],'all')];
colorLimits = [-3e-3, 3e-3];
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_0m(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha 0 %s to %s (UTC)',Time.date_ts(1),Time.date_ts(end)))
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xlim([datenum(Time.date_ts(1)) datenum(Time.date_ts(end))]) 
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
hold off
subplot(3,1,2)
alpha_total_rawm = Alpha.alpha_total_raw .* SNRm .* cloud_SDm_above;
imAlpha=ones(size(alpha_total_rawm(1:ind_km_max,:)));
imAlpha(isnan(alpha_total_rawm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_total_rawm(1:ind_km_max,:),[],'all'), max(alpha_total_rawm(1:ind_km_max,:),[],'all')];
colorLimits = [min(alpha_total_rawm(1:ind_km_max,:),[],'all')/5, max(alpha_total_rawm(1:ind_km_max,:),[],'all')/5];
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),alpha_total_rawm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha raw'))
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
hold off
subplot(3,1,3)
imAlpha=ones(size(Alpha.alpha_totalm(1:ind_km_max,:)));
imAlpha(isnan(Alpha.alpha_totalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(Alpha.alpha_totalm(1:ind_km_max,:),[],'all'), max(Alpha.alpha_totalm(1:ind_km_max,:),[],'all')];
colorLimits = [-2e-4, 10e-4];
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_totalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha total'))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
hold off

%=2D BSR
figure(7473)
BSRm = HSRL.BSR.*cloud_SDm_above.*SNRm;
BSRm = HSRL.BSR;
imagesc(datenum(Time.date_ts),Range.rkm,BSRm)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
hold off
colorbar
colormap(flipud(hot))
set(gca,'ColorScale','log')
caxis([1 50])
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
ylabel('Range (m)')
set(gca, 'YDir','normal')
title('BSR')

%=T sonde and DIAL profiles
figure(685438)
plot(240:300,240:300,'--k')
hold on
plot(242:302,240:300,'k')
plot(238:298,240:300,'k')
for ii = 1:size(Sonde.sonde_ind,2)
    plot(diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,ii))),Sonde.T_sonde(:,ii))
end
hold off
xlabel('T DIAL')
ylabel('T Sonde')

%=T Sonde and DIAL scatter plot
figure(67755)
[~,sondeCut]=min(abs(Range.rm-250));
[~,sondeCut]=min(abs(Range.rm-150));
%[~,sondeCut]=min(abs(Range.rm-1000));
tempProb = zeros(303-278+1);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:length(278:303)
        for j = 1:length(278:303)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(278+i-1)+1)&(diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(278+i-1)-0)& (Sonde.T_sonde(sondeCut:end,ii)<(278+j-1)+1)&(Sonde.T_sonde(sondeCut:end,ii)>=(278+j-1)-0));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc((278:303)-273,(278:303)-273,tempProb,'AlphaData',imAlpha)
colormap(flipud(hot))%colormap hot
%set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
% plot(239:302-273.15,239:302-273.15,'k','LineWidth',2)
% plot(240:303-273.15,237:300-273.15,'k--','LineWidth',2)
% plot(237:300-273.15,240:303-273.15,'k--','LineWidth',2)
plot((278:303)-273,(278:303)-273,'k--')
plot((278:303)-273-2,(278:303)-273,'k')
plot((278:303)-273+2,(278:303)-273,'k')
grid on
hold off
% xlim([239 302]-273.15)
% ylim([239 302]-273.15)
set(gca, 'YDir','normal')
xlabel('Temp retrieval (^oC)')
ylabel('Sonde (^oC)')

%=WV Sonde and DIAL scatter plot
figure(6775551)
[~,sondeCut]=min(abs(Range.rm-250));
[~,sondeCut]=min(abs(Range.rm-150));
minEdge =0;
maxEdge = 12;
binsPer = 4;
maxEdge = maxEdge*binsPer;
tempProb = zeros(maxEdge-minEdge+1);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:length(minEdge:maxEdge)
        for j = 1:length(minEdge:maxEdge)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(AbsHumm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(minEdge+i-1)/binsPer+1)&(diag(AbsHumm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(minEdge+i-1)/binsPer-0)& (Sonde.AbsHum(sondeCut:end,ii)<(minEdge+j-1)/binsPer+1)&(Sonde.AbsHum(sondeCut:end,ii)>=(minEdge+j-1)/binsPer-0));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,tempProb,'AlphaData',imAlpha)
%colormap(flipud(hot))%colormap hot
%set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
plot((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,'k','LineWidth',2)
grid on
hold off
set(gca, 'YDir','normal')
title('WV total correction')
xlabel('WV retrieval (g/m^3)')
ylabel('Sonde (g/m^3)')

%=WV Sonde and DIAL scatter plot
figure(6775552)
[~,sondeCut]=min(abs(Range.rm-250));
[~,sondeCut]=min(abs(Range.rm-150));
minEdge =0;
maxEdge = 12;
binsPer = 4;
maxEdge = maxEdge*binsPer;
tempProb = zeros(maxEdge-minEdge+1);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:length(minEdge:maxEdge)
        for j = 1:length(minEdge:maxEdge)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(AbsHum0m(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(minEdge+i-1)/binsPer+1)&(diag(AbsHum0m(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(minEdge+i-1)/binsPer-0)& (Sonde.AbsHum(sondeCut:end,ii)<(minEdge+j-1)/binsPer+1)&(Sonde.AbsHum(sondeCut:end,ii)>=(minEdge+j-1)/binsPer-0));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,tempProb,'AlphaData',imAlpha)
%colormap(flipud(hot))%colormap hot
%set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
plot((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,'k','LineWidth',2)
grid on
hold off
set(gca, 'YDir','normal')
title('WV 0')
xlabel('WV retrieval (g/m^3)')
ylabel('Sonde (g/m^3)')

figure(67757)
%[~,sondeCut]=min(abs(Range.rm-250));
%[~,sondeCut]=min(abs(Range.rm-1000));
tempProb = zeros(303-278+1);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:length(278:303)
        for j = 1:length(278:303)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(Model.T(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(278+i-1)+1)&(diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(278+i-1)-0)& (Sonde.T_sonde(sondeCut:end,ii)<(278+j-1)+1)&(Sonde.T_sonde(sondeCut:end,ii)>=(278+j-1)-0));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc((278:303)-273,(278:303)-273,tempProb,'AlphaData',imAlpha)
colormap(flipud(hot))%colormap hot
set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
% plot(239:302-273.15,239:302-273.15,'k','LineWidth',2)
% plot(240:303-273.15,237:300-273.15,'k--','LineWidth',2)
% plot(237:300-273.15,240:303-273.15,'k--','LineWidth',2)
plot((278:303)-273,(278:303)-273,'k--')
plot((278:303)-273-2,(278:303)-273,'k')
plot((278:303)-273+2,(278:303)-273,'k')
grid on
hold off
% xlim([239 302]-273.15)
% ylim([239 302]-273.15)
set(gca, 'YDir','normal')
xlabel('lapse rate (K)')
ylabel('Sonde (K)')

%=T Sonde-DIAL histogram
figure(67756)
hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
for ii = 1:size(Sonde.sonde_ind,2)
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Sonde.T_sonde(jj,ii)-Temperature.T_finalm(jj,Sonde.sonde_ind(jj,ii));
        int=int+1;
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('Histogram T_{sonde}-T_{DIAL}\n Mean %f, std %f',meanHist,stdHist))
end
xlabel('\Delta T Sonde-MPD (^oC)')
xlim([-10 10])
ylabel('Occurrences')

figure(67758)
hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;

gg = cloud_SDm_above.*SNRm;
gg(gg<=0)=nan;
ModelTm = Model.T;
ModelTm(isnan(gg)) = NaN;          % Replace mask with NaNs
for ii = 1:size(Sonde.sonde_ind,2)
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Sonde.T_sonde(jj,ii)-ModelTm(jj,Sonde.sonde_ind(jj,ii));
        int=int+1;
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('Histogram T_{sonde}-T_{Lapse}\n Mean %f, std %f',meanHist,stdHist))
end

figure(677564)
hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
for ii = 1:size(Sonde.sonde_ind,2)
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=(Sonde.WV_sonde(jj,ii)-N_wvm(jj,Sonde.sonde_ind(jj,ii)))./Sonde.WV_sonde(jj,ii)*100;
        int=int+1;
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',5)
    title(sprintf('Histogram WV_{sonde}-WV_{DIAL} %%\n Mean %f, std %f',meanHist,stdHist))
end
%xlabel('\Delta T Sonde-MPD (^oC)')
xlim([-40 40])
ylabel('Occurrences')

figure(677565)
hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
for ii = 1:size(Sonde.sonde_ind,2)
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=(Sonde.WV_sonde(jj,ii)-N_wv0m(jj,Sonde.sonde_ind(jj,ii)))./Sonde.WV_sonde(jj,ii)*100;
        int=int+1;
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',5)
    title(sprintf('Histogram WV_{sonde}-WV0_{DIAL}%%\n Mean %f, std %f',meanHist,stdHist))
end
%xlabel('\Delta T Sonde-MPD (^oC)')
xlim([-40 40])
ylabel('Occurrences')

%=FFT of Counts
% figure(7726)
% subplot(2,2,1)
% Fs = 1/Range.rangeBin;
% T = Range.rangeBin;
% L = length(Range.rm);
% f = Fs*(0:(L/2))/L;
% Y = fft(Counts.o2on(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft(Counts.o2off(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(Counts.o2on_mol(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(Counts.o2off_mol(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% title('FFT of counts in range at plot point')
% legend('on','off','on mol','off mol')
% xlabel('f (1/m)')
% ylabel('|P1(f)|')
% %figure(722897)
% subplot(2,2,2)
% Y = fft(fillmissing(Alpha.alpha_total(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft(fillmissing(Alpha.alpha_total_raw(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(fillmissing(Sonde.absorption_sonde{sonde_index}(:,:),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% legend('Alpha total','alpha total raw','alpha sonde')
% title('FFT of absorption in range dimension at plot point')
% xlabel('Freqency (1/m)')
% %figure(722898)
% subplot(2,2,3)
% Y = fft(fillmissing(Temperature.T_final_tests(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft(fillmissing(Temperature.T_final_test(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(fillmissing(Sonde.T_sonde(:,sonde_index),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% legend('T total','T total raw','T sonde')
% title('FFT of Temperature in range dimension at plot point')
% xlabel('Freqency (1/m)')
% %figure(722899)
% subplot(2,2,4)
% Y = fft([ones(1,8) zeros(1,Range.i_range-8)]);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft([ones(1,20) zeros(1,Range.i_range-20)]);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% legend('8 length','10 length')

%=N bins over time
figure(66572)
bar(datenum(Time.date_ts),Counts.NBins)
ylabel('# of bins summed')
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
% xline(738397.94)
% xline(738398.17)
% xline(738400.65)
% xline(738401.18)
% % xline(738407.17)
% % xline(738407.63)
% % xline(738408.93)
% % xline(738409.72)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')

%=Counts and corrections
figure(6)
plot(Range.rm,diag(Counts.o2on(:,p_point)))
hold on
plot(Range.rm,diag(Counts.o2off(:,p_point)))
plot(Range.rm,diag(Counts.o2on_mol(:,p_point)))
plot(Range.rm,diag(Counts.o2off_mol(:,p_point)))
%plot(Range.rm,diag(Counts.delta(:,p_point)),'.-')
%plot(Range.rm,diag(Counts.deltaModel(:,p_point)),'.-')
%plot(Range.rm,diag(Counts.o2on(:,p_point)-Counts.delta(:,p_point)),'--')
%plot(Range.rm,diag(Counts.o2off(:,p_point)-Counts.delta(:,p_point)),'--')
ylim([-200 1000])
%load('AfterPulse.mat')
%plot(Range.rm,pulseON(4:160+3))
%plot(Range.rm,[pulseON(4:160); zeros(3,1)])
% plot(Range.rm,diag(Counts.o2on(:,p_point)-pulseON(4:133+3)),'.-')
% plot(Range.rm,diag(Counts.o2off(:,p_point)-pulseON(4:133+3)),'.-')
%%plot(Range.rm,diag(Counts.o2on(:,p_point)-[pulseON(4:160); zeros(3,1)]),'.-')
%%plot(Range.rm,diag(Counts.o2off(:,p_point)-[pulseON(4:160); zeros(3,1)]),'.-')
% plot(Range.rm,diag(Counts.o2on_mol(:,p_point)))
% plot(Range.rm,diag(Counts.o2off_mol(:,p_point)))
%plot(Range.rm(1:end-1),-diff(Counts.o2off(:,p_point(1))))
%%plot(Range.rm(1:end-1),-diff(Counts.o2off(:,p_point(1)))-[pulseON(4:159); zeros(3,1)])

% plot(Range.rm,Model.N_on(:,p_point(1)))
% plot(Range.rm,Model.N_off(:,p_point(1)))
% plot(Range.rm,Model.N_on_pulse(:,p_point(1)),'--')
% plot(Range.rm,Model.N_off_pulse(:,p_point(1)),'--')
grid on
hold off
%legend('On','off','\Delta','On -\Delta','Off -\Delta','Measured afterpulse')
legend('On','off','On mol','off mol')
title(sprintf(['Counts\n' datestr(Time.date_ts(p_point(1)))]))

% figure(666)
% plot(Range.rm,Model.OverlapOn(:,p_point(1)));
% hold on
% plot(Range.rm,Model.OverlapOff(:,p_point(1)));
% plot(Range.rm,Model.OverlapOn_pulse(:,p_point(1)),'--');
% plot(Range.rm,Model.OverlapOff_pulse(:,p_point(1)),'--');
% hold off
% title('Model counts overlap')
% legend('on','off','on pulse','off pulse')

%=Wavelength overlap correction model
% figure(555)
% plot(Range.rkm,Model.O_on_O_off(:,p_point(1)-10:p_point(1)+10))
% grid on
% title('O_{on}(r)/O_{off}(r) correction')
% xlabel('Range km')

%=Transmission
figure(90)
subplot(2,1,1)
lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(Model.transmission.^2);
lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(Sonde.trasmission_sonde{sonde_index}.^2);
plot(Range.rm,lnOonOoff(:,p_point(1)))
hold on
for time = 1:length(Time.ts)
        % Set fit exclusion zones
        Temperature.exclusion(:,time) = Range.rm<1500 | Range.rm>2500 | SNRm(:,time)==0 | cloud_SDm_above(:,time) ~= 1;
        logicalExc(:,time) = ~logical(Temperature.exclusion(:,time)); 
        V = [Range.rm(logicalExc(:,time)) ones(length(Range.rm(logicalExc(:,time))),1)];
        %calculate polynomial coefficients     
        p = V\Alpha.alpha_0(logicalExc(:,time),time);
        alpha_fit_slope(:,time) = p(1);
        alpha_fit_int(:,time) = p(2);
        alpha_fit(:,time) = alpha_fit_int(:,time) + Range.rm.*alpha_fit_slope(:,time);
end
transmission = exp(-cumtrapz(Range.rm,alpha_fit));
lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(transmission.^2);
%lnOonOoff = lnOonOoff - mean(lnOonOoff(40:67,:),1);
lnOonOoff = lnOonOoff - mean(lnOonOoff(14:27,:),1);
plot(Range.rm,lnOonOoff(:,p_point(1)),'--')
grid on
legend('model lnOonOoff','Fit lnOonOoff')
hold off

%=log on and off
%figure(91)
subplot(2,1,2)
lnOonOoff=exp(lnOonOoff);
plot(Range.rm,lnOonOoff(:,p_point),'--')
%plot(Range.rm,lnOonOoff(:,p_point)./Corr(:,p_point))
title('lnOonOoff')
grid on
hold off

figure(9090)
imagesc(Time.thr,Range.rm/1000,lnOonOoff)
hold on
xline(Time.thr(p_point(1)),'K');
hold off
set(gca, 'YDir','normal')
caxis([.9 1.1])
colorbar
colormap(redblue(64))
title('Counts.delta')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
           
% figure
% plot(Range.rm,log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1))))
% hold on
% plot(Range.rm,log(Model.transmission(:,p_point(1)).^2))
% 
% %TT = Stemp(minJJ) + lapse(minII)/1000 .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
% %PP = Spress(minKK) .* (Stemp(minJJ)./TT).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
% %alpha = absorption_O2_770_model(TT,PP,Spectrum.nu_online(1),Model.WV(:,p_point)); %[m-1] Funcrtion to calculate theoretical absorption
% %Transmission = exp(-cumtrapz(Range.rm,alpha));
% plot(Range.rm,log(Transmission(:,p_point(1)).^2)+const(p_point(1)))
% hold off
% legend('Counts','model','fit')
% 
% figure
% plot(Range.rm,exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Transmission(:,p_point(1)).^2)-const(p_point(1))))
% hold on
% plot(Range.rm,exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Model.transmission(:,p_point(1)).^2)))
% %plot(Range.rm,smooth(exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Transmission(:,p_point(1)).^2)-const(p_point(1))),5))
% plot(Range.rm,overlapCorr(:,p_point))
% hold off
% legend('Fit','model')


% figure(3567)
% plot(diag(Alpha.alpha_total_raw(:,p_point)-Alpha.alpha_0(:,p_point)),Range.rkm)
% xlim([-2e-4 2e-4])
% ylim([0 5])
% grid on
% title(sprintf(['alpha final - alpha 0\n' datestr(Time.date_ts(p_point(1)))]))
% ylabel('Range (km)')
% xlabel('Absorption (m^{-1})')






%=bg
if ishandle(6488)
    figure(6488)
    clf(6488)
end
figure(6488)
line(Time.thr,Counts.bg_o2on-Counts.bg_o2off,'Color','r')
line(Time.thr,Counts.bg_o2on_mol-Counts.bg_o2off_mol,'Color','g')
line(Time.thr,Counts.bg_o2on_mol-Counts.bg_o2on,'Color','g','LineStyle','--')
line(Time.thr,Counts.bg_o2off_mol-Counts.bg_o2off,'Color','r','LineStyle','--')
ylim([-5 5]);
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('Online Background - offline background')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(Time.thr,Temperature.T_finalm(9,:),'Parent',ax2,'Color','k')
line(Time.thr,Temperature.T_final_test(13,:),'Parent',ax2,'Color','b')
ylabel('Temperature horizontal')

if ishandle(64888)
    figure(64888)
    clf(64888)
end
figure(64888)
line(Time.thr,smooth(-Counts.bg_o2on+Counts.bg_o2off,10),'Color','r')
ylim([-5 5])
%ylim([-25 25])
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('-Online Background + offline background')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(Time.thr,smooth(Alpha.alpha_0(4,:),10),'Parent',ax2,'Color','k')
line(Time.thr,smooth(Alpha.alpha_0(9,:),10),'Parent',ax2,'Color','b')
ylim([0 5e-4])
grid on
ylabel('alpha 0')

figure(3927)
plot(Data.Power.O2Online.TimeStamp,Data.Power.O2Online.LaserPower)
hold on
plot(Data.Power.O2Offline.TimeStamp,Data.Power.O2Offline.LaserPower)
%plot(Data.Power.O2Online.TimeStamp,smooth(Data.Power.O2Online.LaserPower,1000))
%plot(Data.Power.O2Offline.TimeStamp,smooth(Data.Power.O2Offline.LaserPower,1000))
xline(Time.ts(p_point(1))/60/60)
hold off
xlim([Time.ts(1) Time.ts(end)]/60/60)
grid on
title('Power Measuremet')
legend('Online','Offline')

figure(337)
subplot(1,2,1)
Model.T_pot = Model.T.*(0.9869232667160128./Model.P).^(0.286);
T_pot = Temperature.T_finalm.*(0.9869232667160128./Model.P).^(0.286);
if ~isempty(Sonde.T_sonde)
    Sonde.T_pot = Sonde.T_sonde.*(0.9869232667160128./Sonde.P_sonde).^(0.286);
end
plot(diag(Model.T_pot(:,p_point)),Range.rkm)
hold on
plot(Sonde.T_pot(:,sonde_index),Range.rkm)
plot(diag(T_pot(:,p_point)),Range.rkm)
hold off
legend('model','measurement','sonde')
grid on
xlabel('Potential Temperature (K)')
ylabel('Range (Range.rkm)')
subplot(1,2,2)
plot(diff(diag(Model.T_pot(:,p_point)))/Range.rangeBin ,Range.rkm(1:end-1))
hold on
plot(diff(diag(T_pot(:,p_point)))/Range.rangeBin,Range.rkm(1:end-1))
plot(diff(Sonde.T_pot(:,sonde_index))/Range.rangeBin,Range.rkm(1:end-1))
xlim([-.15 .15])
hold off
legend('model','measurement','sonde')
grid on
xlabel('Potential Temperature lapse rate (K/m)')
ylabel('Range (Range.rkm)')


figure(338)
r_max_plot = 5; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
r_min_plot = .5; % Max range to plot [km]
[~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
colorLimits = [min(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1
colorLimits =[290 340];
imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),T_pot(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits) %Plot
hold on
%sonde lines
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
%yline(.6)
for ii = 1:size(Sonde.sonde_ind,2)
    plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
end
% % xline(738397.94)
% % xline(738398.17)
% % xline(738400.65)
% % xline(738401.18)
xline(738407.17,'LineWidth',2)
xline(738407.63,'LineWidth',2)
xline(738408.93,'LineWidth',2)
xline(738409.72,'LineWidth',2)
colors = colormap; % Set colors to colormap
colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colormap(colors); % Make colormap new colors
set(gca, 'YDir','normal') % Set y axis increasing
set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
xticks(datenum(Time.date_ts(1)):Format.tickSpacing:datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title(sprintf('Potential Temperature retrieval')) %Add title
xlabel('Time (UTC hours)')
ylabel('Range (km)')
title(colorbar,'Temperature (K)') % Add title to colorbar
hold off

% figure(339)
% r_max_plot = 5; % Max range to plot [km]
% [~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
% r_min_plot = .5; % Max range to plot [km]
% [~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
% imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max-1,:))); % initalize transpernecy matrix
% imAlpha(isnan(Temperature.T_finalm(1:ind_km_max-1,:)))=0; % Set AlphaData to not plot NaNs
% imAlpha(cloud_SDm(1:ind_km_max-1,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
% %imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transoierency so that cloud mask can be seen
% %colorLimits = [min(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1
% colorLimits =[-.5 .5];
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),diff(T_pot(1:ind_km_max,:),1,1),'AlphaData',imAlpha,colorLimits) %Plot
% hold on
% %sonde lines
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% %yline(.6)
% for ii = 1:size(Sonde.sonde_ind,2)
%     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% end
% % % xline(738397.94)
% % % xline(738398.17)
% % % xline(738400.65)
% % % xline(738401.18)
% xline(738407.17,'LineWidth',2)
% xline(738407.63,'LineWidth',2)
% xline(738408.93,'LineWidth',2)
% xline(738409.72,'LineWidth',2)
% colors = colormap; % Set colors to colormap
% colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
% colormap(colors); % Make colormap new colors
% set(gca, 'YDir','normal') % Set y axis increasing
% set(gca,'color',[1 1 1]);% Color background white
% colorbar % Add colorbar
% xticks(datenum(Time.date_ts(1)):Format.tickSpacing:datenum(Time.date_ts(end)))
% xtickangle(Format.tickAngle)
% datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
% title(sprintf('Potential Temperature retrieval')) %Add title
% xlabel('Time (UTC hours)')
% ylabel('Range (km)')
% title(colorbar,'Temperature (K)') % Add title to colorbar
% hold off

figure(2828)
title('time')
subplot(1,4,1)
semilogx(Counts.Poissonthin.timeWidthon*Time.t_step,Range.rkm,'.')
xlim([1 10^4])
grid on
title('on')
ylabel('Range')
subplot(1,4,2)
semilogx(Counts.Poissonthin.timeWidthoff*Time.t_step,Range.rkm,'.')
xlim([1 10^4])
grid on
title('off')
subplot(1,4,3)
semilogx(Counts.Poissonthin.timeWidthon_mol*Time.t_step,Range.rkm,'.')
xlim([1 10^4])
grid on
title('on mol')
subplot(1,4,4)
semilogx(Counts.Poissonthin.timeWidthoff_mol*Time.t_step,Range.rkm,'.')
xlim([1 10^4])
grid on
title('of mol')

figure(2829)
title('time')
subplot(4,1,1)
semilogy(Time.thr,Counts.Poissonthin.rangeWidthon*Range.rangeBin,'.')
ylim([1 300])
xlim([Time.thr(1) Time.thr(end)])
grid on
title('on')
subplot(4,1,2)
semilogy(Time.thr,Counts.Poissonthin.rangeWidthoff*Range.rangeBin,'.')
ylim([1 300])
xlim([Time.thr(1) Time.thr(end)])
grid on
title('off')
subplot(4,1,3)
semilogy(Time.thr,Counts.Poissonthin.rangeWidthon_mol*Range.rangeBin,'.')
ylim([1 300])
xlim([Time.thr(1) Time.thr(end)])
grid on
title('on_mol')
subplot(4,1,4)
semilogy(Time.thr,Counts.Poissonthin.rangeWidthoff_mol*Range.rangeBin,'.')
ylim([1 300])
xlim([Time.thr(1) Time.thr(end)])
grid on
title('of_mol')
xlabel('Time UTC')

