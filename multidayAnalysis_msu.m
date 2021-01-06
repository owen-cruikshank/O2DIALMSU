%%D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\6_24_20\20190417.mat'

%load('D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\7_20_20\20200717.mat',...
%'alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','o2on','o2off','o2on_mol','o2off_mol'...
%    ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
%    ,'absorption_f','InsideCell','OutsideCell','TSOA','RoomTemp','OnlineWaveDiff','OfflineWaveDiff'...
%    ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol')

load('C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\7_20_20\20200717.mat',...
'alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','o2on','o2off','o2on_mol','o2off_mol'...
    ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
    ,'absorption_f','InsideCell','OutsideCell','TSOA','RoomTemp','OnlineWaveDiff','OfflineWaveDiff'...
    ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol')
    
alpha_totalm_O = alpha_totalm;
alpha_0_O = alpha_0;
T_O = T;
T_finalm_O = T_finalm;
BSR_O = BSR;
ts_O=ts;
cloud_SDm_O=cloud_SDm;
rm_O=rm;
%rkm_O=rkm;
o2on_O=o2on;
o2off_O=o2off;
o2on_mol_O=o2on_mol;
o2off_mol_O=o2off_mol;
Ts_fit_O=Ts_fit;
L_fit_sm_test_O=L_fit_sm_test;
bg_o2off_O=bg_o2off;
bg_o2on_O=bg_o2on;
bg_o2off_mol_O=bg_o2off_mol;
bg_o2on_mol_O=bg_o2on_mol;
alpha_total_raw_O=alpha_total_raw;
absorption_O=absorption;
absorption_f_O=absorption_f;
%mask_nodata_O=mask_nodata;
o2on_intp_O=o2on_intp;
o2off_intp_O=o2off_intp;
o2on_intp_mol_O=o2on_intp_mol;
o2off_intp_mol_O=o2off_intp_mol;
%o2on_noise_O=o2on_noise;
o%2on_noise_mol_O=o2on_noise_mol;
%o2off_noise_O=o2off_noise;
%o2off_noise_mol_O=o2off_noise_mol;
%date_ts_N_O=date_ts_N;
%data_col_real_O=data_col_real;
%date_ts_O=date_ts;
%WV_O=WV;
%Ts_O=Ts;
%exclusion_O=exclusion;
%T_final_test_O=T_final_test;
%T_real_O=T_real;
%absorption_sonde_O=absorption_sonde;

    %%%D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\6_24_20\20190420.mat'
load('C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\7_20_20\20200720.mat','alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','o2on','o2off','o2on_mol','o2off_mol'...
    ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
    ,'absorption_f','InsideCell','OutsideCell','TSOA','RoomTemp','OnlineWaveDiff','OfflineWaveDiff'...
    ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol')
    
alpha_totalm = [alpha_totalm_O alpha_totalm];
 alpha_0=[alpha_0_O  alpha_0];
T=[T_O  T];
T_finalm=[T_finalm_O  T_finalm];
BSR=[BSR_O  BSR];
cloud_SDm=[cloud_SDm_O cloud_SDm];
%rm=[rm_O rm];
%rkm =[rkm_O rkm];
o2on=[o2on_O o2on];
o2off=[o2off_O o2off];
o2on_mol=[o2on_mol_O o2on_mol];
o2off_mol=[o2off_mol_O o2off_mol];
Ts_fit=[Ts_fit_O Ts_fit];
L_fit_sm_test=[L_fit_sm_test_O L_fit_sm_test];
bg_o2off=[bg_o2off_O bg_o2off];
bg_o2on=[bg_o2on_O bg_o2on];
bg_o2off_mol=[bg_o2off_mol_O bg_o2off_mol];
bg_o2on_mol=[bg_o2on_mol_O bg_o2on_mol];
alpha_total_raw=[alpha_total_raw_O alpha_total_raw];
absorption=[absorption_O absorption];
absorption_f=[absorption_f_O absorption_f];
mask_nodata=[mask_nodata_O mask_nodata];
o2on_intp=[o2on_intp_O o2on_intp];
o2off_intp=[o2off_intp_O o2off_intp];
o2on_intp_mol=[o2on_intp_mol_O o2on_intp_mol];
o2off_intp_mol=[o2off_intp_mol_O o2off_intp_mol];
o2on_noise=[o2on_noise_O o2on_noise];
o2on_noise_mol=[o2on_noise_mol_O o2on_noise_mol];
o2off_noise=[o2off_noise_O o2off_noise];
o2off_noise_mol=[o2off_noise_mol_O o2off_noise_mol];
WV=[WV_O WV];
Ts=[Ts_O Ts];
exclusion=[exclusion_O exclusion];
T_final_test=[T_final_test_O T_final_test];
T_real=[T_real_O T_real];
absorption_sonde=[absorption_sonde_O absorption_sonde];


date_ts_N=[date_ts_N_O date_ts_N];
date_ts=[date_ts_O date_ts];

ts=[ts_O (ts+ts_O(end))];
data_col_real=[data_col_real_O (data_col_real+length(ts_O))];

%remove nans
sondeNans = isnan(data_col_real);
data_col_real = data_col_real(~sondeNans);
T_real = T_real(:,~sondeNans);



%%

%===============
%=== Figures ===
%===============


%Tick spacing
tickHours = 8;
tickMin = 0;
date2=datetime(2019,4,20,tickHours,tickMin,0);
date1=datetime(2019,4,20,0,0,0);
tickSpacing = datenum(date2)-datenum(date1);
tickAngle = 45;



t_index = 21;    % Choose which measurement to compare. Max value here is number of "real" measurements
p_point = data_col_real(t_index);



figure(728590)
plot(datenum(date_ts_N),permute(L_fit_sm_test(1,:,:),[3 2 1]))
xline(datenum(date_ts_N(p_point)));
grid on

xlim([datenum(date_ts_N(1)) datenum(date_ts_N(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
xtickangle(tickAngle)
ylim([-11e-3 -3e-3]) 
title('Fitted lapse rate')
xlabel('Time (UTC)')
ylabel('Lapse rate (K/km)')

figure(728591)
%plot(ts/60/60,Ts - permute(Ts_fit(1,:,:),[3 2 1]))
plot(ts/60/60,T(1,:) - permute(Ts_fit(1,:,:),[3 2 1]))
xline(ts(p_point)/60/60);
grid on
title('Surface temp - fitted surface temp')
xlabel('time (hr)')
ylabel('\deltaT (K)')




figure(99492)
plot(BSR(:,p_point),rkm)
legend('BSR')
title({'Backscatter Ratio';datestr(date_ts_N(p_point))})
ylabel('Range (km)')





figure(99499)
plot(WV(:,p_point),rkm)
xlabel('WV (molec/m^2')
ylabel('Range (km)')
title({'Water vapor profile';datestr(date_ts_N(p_point))})







figure(884)
plot(T(:,p_point),rkm)
hold on
%plot(T(:,p_point)+2,rkm)
%plot(T(:,p_point)-2,rkm)
%plot(Ts(p_point),0,'+')
%plot(Ts_fit(1,p_point,end),0,'+')
plot(T_finalm(:,p_point),rkm)
%plot(T_final(:,p_point),rkm)
plot(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end),rkm,'--')
%plot(T_final(:,p_point),rm)
exclusion(exclusion==0)=NaN;
plot(T_final_test(:,p_point).*exclusion(:,p_point,end),rkm,'*')
plot(T_real(:,t_index),rkm,'.-')
plot(T_real(:,t_index)-2,rkm,'.-')
plot(T_real(:,t_index)+2,rkm,'.-')
hold off
grid on
xlabel('Temperature (K)')
ylabel('Range (km)')
ylim([-inf 4])
title({'Temperature profile';datestr(date_ts_N(p_point))})
legend('Temp guess','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde','Sonde-2K','Sonde+2K''Location','southwest')

% Single time temperature deviation
figure(885)
plot(T_real(:,t_index)-T_finalm(:,p_point),rkm)
%title(sprintf('Temperature Deviation  %s (UTC)',date_time))
line([0 0],[0 6],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--')
line([-1 -1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.')
line([1 1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
line([2 2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
title({'\Delta T (T_{sonde} - T_{retrieved}) (K)';datestr(date_ts_N(p_point))})
ylabel('Range (km)')
legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','Southwest')
ylim([0 4])
xlim([-10 10])




figure(9985)
bg_o2offm = bg_o2off.*mask_nodata(1,:);
bg_o2offm(bg_o2offm <= 0) = NaN;
[bg_o2offm,TF]=rmmissing(bg_o2offm,2);
bg_o2onm = bg_o2on.*mask_nodata(1,:);
bg_o2onm(bg_o2onm <= 0) = NaN;
bg_o2onm=rmmissing(bg_o2onm,2);
bg_o2off_molm = bg_o2off_mol.*mask_nodata(1,:);
bg_o2off_molm(bg_o2off_molm <= 0) = NaN;
bg_o2off_molm=rmmissing(bg_o2off_molm,2);
bg_o2on_molm = bg_o2on_mol.*mask_nodata(1,:);
bg_o2on_molm(bg_o2on_molm <= 0) = NaN;
bg_o2on_molm=rmmissing(bg_o2on_molm,2);

semilogy(datenum(date_ts_N(~TF)),bg_o2offm)
hold on
semilogy(datenum(date_ts_N(~TF)),bg_o2onm)
semilogy(datenum(date_ts_N(~TF)),bg_o2off_molm)
semilogy(datenum(date_ts_N(~TF)),bg_o2on_molm)
legend('Offline','Online','Offline molecular','Online molecular')
xlim([datenum(date_ts_N(1)) datenum(date_ts_N(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
xtickangle(tickAngle)
title('Background')
xlabel('Time UTC')
ylabel('photons')
hold off
grid on






figure(886)
% if length(span_days)==1
%     title_temp = append('Temperature profile for ' , datestr(span_days));
% else
%     title_temp = 'Temperature profile for ';
%     for i = 1:length(span_days)
%         title_temp = append(title_temp , datestr(span_days(i)));
%     end
% end
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(T_finalm(1:ind_km_max,:)));
imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(T_finalm(1:ind_km_max,:)-1,[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:length(data_col_real)
    if ~isnan(data_col_real(i))
        xline(datenum(date_ts_N(data_col_real(i))));
    end
end
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('Temperature %s to %s (UTC)',date_ts(1),date_ts(end)))
%title(sprintf('Temperature retrieval\n%s(UTC)',span_days(1)))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'Temperature (K)')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off
%title('Temperature retrieval')
%title(file)


figure(98872)
subplot(2,1,1)
imagesc(datenum(date_ts),rkm(1:ind_km_max),o2on(1:ind_km_max,:))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
set(gca,'colorscale','log')
set(gca, 'YDir','normal')
subplot(2,1,2)
imagesc(datenum(date_ts),rkm(1:ind_km_max),o2off(1:ind_km_max,:))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
set(gca,'colorscale','log')
set(gca, 'YDir','normal')

figure(98873)
imagesc(datenum(date_ts),rkm(1:ind_km_max),BSR(1:ind_km_max,:))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
set(gca,'colorscale','log')
set(gca, 'YDir','normal')




figure(4)
plot(alpha_0(:,p_point),rkm)
hold on
plot(absorption(:,p_point),rkm)
plot(alpha_totalm(:,p_point),rkm)
plot(alpha_total_raw(:,p_point),rkm)
plot(absorption_sonde(:,t_index),rkm,'-.')
legend('Measured zeroth order absorption','Theoretical absorption from Tg','Pertabative absorption','Alpha 0+1+2','Sonde')%,'theoretical absorption online-offine')
hold off
grid on
ylim([-inf 4])
title({'O2 absorption profile';datestr(date_ts_N(p_point))})
ylabel('Range (km)')
xlabel('Absorption m^{-1}')





figure(8)
plot(rm,o2on(:,p_point).*rm.^2)
hold on
plot(rm,o2off(:,p_point).*rm.^2)
%title('Averaged photon returns multiplied by r^2')
title({'Averaged counts returns multiplied by r^2';datestr(date_ts_N(p_point))})
legend('Online','Offline')
xlabel('Range m')
ylabel('Counts * r^2')
grid on
hold off

figure(7)
i_range = length(rm);
plot(rm,o2on(:,p_point))
hold on
plot(rm,o2off(:,p_point))
plot(rm,o2on_mol(:,p_point))
plot(rm,o2off_mol(:,p_point))

plot(rm,o2on_noise(1:i_range,p_point),'--')
plot(rm,o2off_noise(1:i_range,p_point),'--')
plot(rm,o2on_noise_mol(1:i_range,p_point),'--')
plot(rm,o2off_noise_mol(1:i_range,p_point),'--')
title({'Averaged counts';datestr(date_ts_N(p_point))})
legend('Online','Offline','Online molecular','Offline molecular')
hold off
grid on
xlabel('Range [m]')
ylabel('Counts')

figure(8)
overlapMinIndex = 20;
%overlap = (o2on./max(o2on(8:end,:),[],1))./(o2on_mol./max(o2on_mol(8:end,:),[],1));
overlap = o2on./o2on_mol;
overlap = overlap./max(overlap(overlapMinIndex:round(end/2),:),[],1);
plot(rm,overlap(:,p_point))

% overlapFile = overlap(:,p_point);
% rmOverlap = rm;
% save('overlapFile4_20_19.mat','overlapFile','rmOverlap');



% % % % figure(9)
% % % % %plot(rm,o2on(:,p_point))
% % % % hold on
% % % % %plot(rm,o2off(:,p_point))
% % % % %plot(rm,o2on_mol(:,p_point))
% % % % %plot(rm,o2off_mol(:,p_point))
% % % % 
% % % % plot(rm_raw,O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data(:,p_point),'--')
% % % % plot(rm_raw,O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data(:,p_point),'--')
% % % % plot(rm_raw,o2off_raw_corrected(:,p_point),'--')
% % % % plot(rm_raw,o2on_raw_corrected(:,p_point),'--')
% % % % 
% % % % ylim([0 500])
% % % % title({'Raw counts';datestr(date_ts_N(p_point))})
% % % % legend('Online','Offline','Online molecular','Offline molecular')
% % % % hold off
% % % % grid on
% % % % xlabel('Range [m]')
% % % % ylabel('Photons')


figure(6773)
%temp diff sonde
hold on
index=1;
for i = 1:length(data_col_real)
    if ~isnan(data_col_real(i))
        plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*')
        fullTFinal(:,index) = T_finalm(:,data_col_real(i));
        T_real_final(:,index) = T_real(:,i);
        index=index+1;
    end

end
plot(270:300,270:300)
plot(272:302,270:300)
plot(268:298,270:300)
xlabel('Temperature Retrieval (K)')
ylabel('Sonde (K)')
title('All sonde measurements')
hold off
grid on

figure(6774)
for i = 1:length(260:305)
    for j = 1:length(260:305)
        
        tempProb(i,j) = nnz((fullTFinal<(260+i-1)+0.5)&(fullTFinal>=(260+i-1)-0.5)& (T_real_final<(260+j-1)+0.5)&(T_real_final>=(260+j-1)-0.5));
    end
end
imagesc(260:305,260:305,tempProb')
colorbar
set(gca, 'YDir','normal')
xlabel('Temp retrieval (K)')
ylabel('Sonde (K)')



figure(6775)
cloud_index = [6 11 12 13 14];
cloud_index = [1 3 4 5 7 8 9];
cloud_index = [1 3 4 5 7 8 9 18 19 23 24 25 26];%4_17to4_22
%cloud_index = [1 2 3 4 10 12 16 17 18 19 20 21 22 23 24 34 35];%5_17 to 5_22
%temp diff sonde
hold on
index=1;
for i = cloud_index
    if ~isnan(data_col_real(i))
        plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*')
        fullTFinal(:,index) = T_finalm(:,data_col_real(i));
        T_real_final(:,index) = T_real(:,i);
        index=index+1;
    end

end
plot(270:300,270:300)
plot(272:302,270:300)
plot(268:298,270:300)
xlabel('Temperature Retrieval (K)')
ylabel('Sonde (K)')
title('Sonde with clouds')
hold off
grid on


figure(6776)
no_cloud_index = [2 6 10 11 12 13 14 15 16 17 20 21 22];%4_17to4_22
%no_cloud_index = [5 6 7 8 9 11 13 14 15 25 26 27 28 29 30 31 32 33 36 37 38 39];%5_17 to 5_22

%temp diff sonde
hold on
index=1;
for i = no_cloud_index
    if ~isnan(data_col_real(i))
        plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*')
        fullTFinal(:,index) = T_finalm(:,data_col_real(i));
        T_real_final(:,index) = T_real(:,i);
        index=index+1;
    end

end
plot(270:300,270:300)
plot(272:302,270:300)
plot(268:298,270:300)
xlabel('Temperature Retrieval (K)')
ylabel('Sonde (K)')
title('Sonde without clouds')
hold off
grid on
