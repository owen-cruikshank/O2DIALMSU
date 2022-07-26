clear all
file='C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\2_18_21\';
%file='D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\2_18_21\';

file='C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\3_29_21\newOwenProgram\';
file='C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\3_29_21\newOwenProgramNoDelta\';
file='C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\3_29_21\newnewOwenProgramNoDelta\';
file='C:\Users\d98b386\OneDrive - Montana State University\research s19\Reports\3_29_21\newnewnewOwenProgramNoDelta\';



files = [20190417 20190418 20190419 20190420 20190421 20190422];
%files = [20190418 20190419 20190420 20190421 20190422];

files = [20200829 20200830 20200831 20200901 20200902 20200903 20200904 20200905 20200906 20200907 20200908];

%file1='D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\6_24_20\20190510.mat';
%file2='D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\6_24_20\20190513.mat';


% load([file num2str(files(1)) '.mat'],'alpha_totalm'...
%         ,'alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
%         ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
%         ,'absorption_f','mask_nodata'...
%         ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
%         ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
%         ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')
    
    load([file num2str(files(1)) '.mat'],'alpha_totalm'...
        ,'alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
        ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
        ,'mask_nodata'...
        ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
        ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
        ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')
    
alpha_totalm_O = alpha_totalm;
alpha_0_O = alpha_0;
T_O = T;
T_finalm_O = T_finalm;
BSR_O = BSR;
ts_O=ts;
cloud_SDm_O=cloud_SDm;
rm_O=rm;
rkm_O=rkm;
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
%absorption_f_O=absorption_f;
mask_nodata_O=mask_nodata;
o2on_intp_O=o2on_intp;
o2off_intp_O=o2off_intp;
o2on_intp_mol_O=o2on_intp_mol;
o2off_intp_mol_O=o2off_intp_mol;
o2on_noise_O=o2on_noise;
o2on_noise_mol_O=o2on_noise_mol;
o2off_noise_O=o2off_noise;
o2off_noise_mol_O=o2off_noise_mol;
date_ts_N_O=date_ts_N;
data_col_real_O=data_col_real;
date_ts_O=date_ts;
WV_O=WV;
Ts_O=Ts;
exclusion_O=exclusion;
T_final_test_O=T_final_test;
T_real_O=T_real;
absorption_sonde_O=absorption_sonde;
    
for j=2:length(files)
% load([file num2str(files(j)) '.mat'],'alpha_totalm'...
%         ,'alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
%         ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
%         ,'absorption_f','mask_nodata'...
%         ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
%         ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
%         ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')
    
    load([file num2str(files(j)) '.mat'],'alpha_totalm'...
        ,'alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
        ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
        ,'mask_nodata'...
        ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
        ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
        ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')
    
alpha_totalm_O = [alpha_totalm_O alpha_totalm];
 alpha_0_O=[alpha_0_O  alpha_0];
T_O=[T_O  T];
T_finalm_O=[T_finalm_O  T_finalm];
BSR_O=[BSR_O  BSR];
cloud_SDm_O=[cloud_SDm_O cloud_SDm];
%rm=[rm_O rm];
%rkm =[rkm_O rkm];
o2on_O=[o2on_O o2on];
o2off_O=[o2off_O o2off];
o2on_mol_O=[o2on_mol_O o2on_mol];
o2off_mol_O=[o2off_mol_O o2off_mol];
Ts_fit_O=[Ts_fit_O Ts_fit];
L_fit_sm_test_O=[L_fit_sm_test_O L_fit_sm_test];
bg_o2off_O=[bg_o2off_O bg_o2off];
bg_o2on_O=[bg_o2on_O bg_o2on];
bg_o2off_mol_O=[bg_o2off_mol_O bg_o2off_mol];
bg_o2on_mol_O=[bg_o2on_mol_O bg_o2on_mol];
alpha_total_raw_O=[alpha_total_raw_O alpha_total_raw];
absorption_O=[absorption_O absorption];
%absorption_f_O=[absorption_f_O absorption_f];
mask_nodata_O=[mask_nodata_O mask_nodata];
o2on_intp_O=[o2on_intp_O o2on_intp];
o2off_intp_O=[o2off_intp_O o2off_intp];
o2on_intp_mol_O=[o2on_intp_mol_O o2on_intp_mol];
o2off_intp_mol_O=[o2off_intp_mol_O o2off_intp_mol];
o2on_noise_O=[o2on_noise_O o2on_noise];
o2on_noise_mol_O=[o2on_noise_mol_O o2on_noise_mol];
o2off_noise_O=[o2off_noise_O o2off_noise];
o2off_noise_mol_O=[o2off_noise_mol_O o2off_noise_mol];
WV_O=[WV_O WV];
Ts_O=[Ts_O Ts];
exclusion_O=[exclusion_O exclusion];
T_final_test_O=[T_final_test_O T_final_test];
T_real_O=[T_real_O T_real];
absorption_sonde_O=[absorption_sonde_O absorption_sonde];



date_ts_N_O=[date_ts_N_O date_ts_N];
date_ts_O=[date_ts_O date_ts];

ts_O=[ts_O (ts+ts_O(end))];
data_col_real_O=[data_col_real_O (data_col_real+length(ts_O)-length(ts))];
end



%%
%remove nans
% % % sondeNans = isnan(data_col_real_O);
% % % for i = 1:length(sondeNans(1,:))
% % % data_col_real_O(:,i) = data_col_real_O(~sondeNans(:,i));
% % % T_real_O(:,i) = T_real_O(:,~sondeNans(1,i));
% % % end


%%
no_cloud_index = [2 6 10 11 12 13 14 15 16 17 20 21 22];%4_17 to 4_22
cloud_index = [1 4 8 9 18 24 25 26];%4_17 to 4_22
all_index = [no_cloud_index cloud_index];
garbage_index = [3 5 7 19 23];

no_cloud_index = [];
cloud_index = [];
all_index = [];
garbage_index = [];


% % no_cloud_index = [ 6 10 11 12 13 14 15 16 17 20 21 22];%4_17 to 4_22
% % cloud_index = [ 8 9 18 24 25 26];%4_17 to 4_22
% % all_index = [no_cloud_index cloud_index];
% % garbage_index = [ 7 19 23];
% % 
% % all_index = all_index-5;


%%
figure(6773)
%temp diff sonde
hold on
index=1;
for i = 1:size(data_col_real_O,2)
    if ~isnan(data_col_real_O(:,i))
        plot(diag(T_finalm_O(:,data_col_real_O(:,i))),T_real_O(:,i),'*')
        fullTFinal_O(:,index) = diag(T_finalm_O(:,data_col_real_O(:,i)));
        T_real_final_O(:,index) = T_real_O(:,i);
        index=index+1;
    end

end
plot(270:300,270:300)
plot(272:302,270:300)
plot(268:298,270:300)
xlabel('Temperature Retrieval (K)')
ylabel('Sonde (K)')
hold off
grid on

% % figure(67755)
% % index = cloud_index;
% % index = all_index;
% % %index = garbage_index;
% % for i = 1:length(260:330)
% %     for j = 1:length(260:330)
% %         
% %         %tempProb(i,j) = nnz((fullTFinal_O(:,index)<(260+i-1)+0.5)&(fullTFinal_O(:,index)>=(260+i-1)-0.5)& (T_real_final_O(:,index)<(260+j-1)+0.5)&(T_real_final_O(:,index)>=(260+j-1)-0.5));
% %         tempProb(i,j) = nnz((fullTFinal_O(:,index)<(260+i-1)+1)&(fullTFinal_O(:,index)>=(260+i-1)-0)& (T_real_final_O(:,index)<(260+j-1)+1)&(T_real_final_O(:,index)>=(260+j-1)-0));
% %     end
% % end
% % tempProb(tempProb==0)=nan;
% % imAlpha=ones(size(tempProb));
% % imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% % imagesc(260:330,260:330,tempProb,'AlphaData',imAlpha)
% % %colormap hot
% % colormap(flipud(hot))
% % set(gca,'ColorScale','log')
% % set(gca,'Color','#D3D3D3')
% % a = colorbar;
% % a.Label.String = 'Occurances';
% % hold on
% % plot(259:302,259:302,'k','LineWidth',2)
% % plot(260:303,257:300,'k--','LineWidth',2)
% % plot(257:300,260:303,'k--','LineWidth',2)
% % hold off
% % xlim([259 302])
% % ylim([259 302])
% % set(gca, 'YDir','normal')
% % xlabel('Temp retrieval (K)')
% % ylabel('Sonde (K)')


% % % figure(3333)
% % % index = all_index;
% % % sondeDiff=reshape(T_real_final_O(:,index)-fullTFinal_O(:,index),numel(T_real_final_O(:,index)),[]);
% % % sondeDiff(abs(sondeDiff)>10)=nan;
% % % sondeDiff = sondeDiff(~isnan(sondeDiff));
% % % histogram(sondeDiff)
% % % 
% % % standardDeviation=std(sondeDiff,'omitnan');
% % % meanHist = mean(sondeDiff,'omitnan');
% % % title('Histogram of (T_{sonde} - T_{DIAL})')
% % % xlabel('(K)')
% % % ylabel('occurcances')
% % % xlim([-9 9])


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



t_index =13;    % Choose which measurement to compare. Max value here is number of "real" measurements
%p_point = data_col_real_O(:,t_index);

time = 184;
[~,p_point] = min(abs(time-ts_O/60/60));

p_point = p_point*ones(size(rm,1),1);

figure(728590)
plot(datenum(date_ts_N_O),permute(L_fit_sm_test_O(1,:,:),[3 2 1]))
%xline(datenum(date_ts_N_O(p_point)));

hold on
plot(datenum(date_ts_N_O([p_point(1) p_point(end)])),[-10 -4]*10^-3)
hold off
grid on

xlim([datenum(date_ts_N_O(1)) datenum(date_ts_N_O(end))]) 
xticks(datenum(date_ts_O(1)): tickSpacing :datenum(date_ts_O(end)))
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
xtickangle(tickAngle)
ylim([-11e-3 -3e-3]) 
title('Fitted lapse rate')
xlabel('Time (UTC)')
ylabel('Lapse rate (K/km)')

figure(728591)
%plot(ts/60/60,Ts - permute(Ts_fit(1,:,:),[3 2 1]))
plot(ts_O/60/60,T_O(1,:) - permute(Ts_fit_O(1,:,:),[3 2 1]))
%xline(ts_O(p_point)/60/60);
hold on
plot(ts_O([p_point(1) p_point(end)])/60/60,[-10 10])
grid on
title('Surface temp - fitted surface temp')
xlabel('time (hr)')
ylabel('\deltaT (K)')




figure(99492)
plot(diag(BSR_O(:,p_point)),rkm_O,'linewidth',2)
legend('BSR')
title({'Backscatter Ratio';datestr(date_ts_N_O(p_point(1)))})
ylabel('Range (km)')





figure(99499)
plot(diag(WV_O(:,p_point)),rkm_O,'linewidth',2)
xlabel('WV (molec/m^2')
ylabel('Range (km)')
title({'Water vapor profile';datestr(date_ts_N_O(p_point(1)))})







figure(884)
plot(diag(T_O(:,p_point)),rkm_O)

%plot(T(:,p_point)+2,rkm)
%plot(T(:,p_point)-2,rkm)
%plot(Ts(p_point),0,'+')
%plot(Ts_fit(1,p_point,end),0,'+')

%hold on
%plot(T_final(:,p_point),rkm)
%plot(diag(L_fit_sm_test_O(:,p_point,end) .* rm_O + Ts_fit_O(1,p_point,end)),rkm_O,'--')
%plot(T_final(:,p_point),rm)
exclusion_O(exclusion_O==0)=NaN;
%plot(diag(T_final_test_O(:,p_point).*exclusion_O(:,p_point,end)),rkm_O,'*')
%%plot(T_real_O(:,t_index),rkm_O,'k-','Linewidth',2)
hold on
%%plot(T_real_O(:,t_index)-2,rkm_O,'b-.','Linewidth',2)
%%plot(T_real_O(:,t_index)+2,rkm_O,'b-.','Linewidth',2)
plot(diag(T_finalm_O(:,p_point)),rkm_O,'r-','Linewidth',2)
hold off
grid on
xlabel('Temperature (K)')
ylabel('Range (km)')
ylim([-inf 5])
title({'Temperature profile';datestr(date_ts_N_O(p_point(1)))})
%legend('Temp guess','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde','Sonde-2K','Sonde+2K','Location','southwest')
%legend('Retrieved Temperature','Sonde','Sonde-2K','Sonde+2K','Location','southwest')
legend('Retrieved Temperature','Temperature guess','Location','southwest')


figure(884)
plot(diag(T_O(:,p_point)),rkm_O,'linewidth',2)

%plot(T(:,p_point)+2,rkm)
%plot(T(:,p_point)-2,rkm)
%plot(Ts(p_point),0,'+')
%plot(Ts_fit(1,p_point,end),0,'+')

%hold on
%plot(T_final(:,p_point),rkm)
%plot(diag(L_fit_sm_test_O(:,p_point,end) .* rm_O + Ts_fit_O(1,p_point,end)),rkm_O,'--','linewidth',2)
%plot(T_final(:,p_point),rm)
exclusion_O(exclusion_O==0)=NaN;
%plot(diag(T_final_test_O(:,p_point).*exclusion_O(:,p_point,end)),rkm_O,'*')
%plot(T_real_O(:,t_index),rkm_O,'k-','Linewidth',2)
hold on
%plot(T_real_O(:,t_index)-2,rkm_O,'b-.','Linewidth',2)
%plot(T_real_O(:,t_index)+2,rkm_O,'b-.','Linewidth',2)
plot(diag(T_finalm_O(:,p_point)),rkm_O,'r-','Linewidth',2)
hold off
grid on
xlabel('Temperature (K)')
ylabel('Range (km)')
ylim([-inf 5])
title({'Temperature profile';datestr(date_ts_N_O(p_point(1)))})
%legend('Temp guess','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde','Sonde-2K','Sonde+2K','Location','southwest')
%legend('Retrieved Temperature','Sonde','Sonde-2K','Sonde+2K','Location','southwest')
%legend('Fitted Lapse Rate','Retrieved Temperature')
legend('Retrieved Temperature','Temperature guess','Location','southwest')

% figure(884)
% %plot(diag(T_O(:,p_point)),rkm_O)
% 
% %plot(T(:,p_point)+2,rkm)
% %plot(T(:,p_point)-2,rkm)
% %plot(Ts(p_point),0,'+')
% %plot(Ts_fit(1,p_point,end),0,'+')
% 
% %hold on
% %plot(T_final(:,p_point),rkm)
% %plot(diag(L_fit_sm_test_O(:,p_point,end) .* rm_O + Ts_fit_O(1,p_point,end)),rkm_O,'--')
% %plot(T_final(:,p_point),rm)
% %exclusion_O(exclusion_O==0)=NaN;
% %plot(diag(T_final_test_O(:,p_point).*exclusion_O(:,p_point,end)),rkm_O,'*')
% %plot(T_real_O(:,t_index),rkm_O,'k-','Linewidth',2)
% %hold on
% %plot(T_real_O(:,t_index)-2,rkm_O,'b-.','Linewidth',2)
% %plot(T_real_O(:,t_index)+2,rkm_O,'b-.','Linewidth',2)
% plot(diag(T_finalm_O(:,p_point))-273,rkm_O,'r-','Linewidth',2)
% %hold off
% grid on
% xlabel('Temperature (C)')
% ylabel('Range (km)')
% ylim([-inf 4])
% title({'Temperature profile';datestr(date_ts_N_O(p_point(1)))})
% %legend('Temp guess','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde','Sonde-2K','Sonde+2K','Location','southwest')
% %legend('Retrieved Temperature','Sonde','Sonde-2K','Sonde+2K','Location','southwest')

figure(8859)
%plot(T_real_O(:,1)-diag(T_finalm_O(:,data_col_real_O(:,1))),rkm_O)
hold on
%for i=2:length(data_col_real_O(1,:))
for i=all_index
    plot(T_real_O(:,i)-diag(T_finalm_O(:,data_col_real_O(:,i))),rkm_O)
end
%title(sprintf('Temperature Deviation  %s (UTC)',date_time))
line([0 0],[0 6],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--')
line([-1 -1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.')
line([1 1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
line([2 2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
title({'\Delta T (T_{sonde} - T_{retrieved}) (K)';datestr(date_ts_N_O(p_point(1)))})
ylabel('Range (km)')
legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','Southeast')
ylim([0 5])
xlim([-15 15])

figure(9949239)
BSRm = BSR_O;
BSRm(~isnan(T_finalm))=nan;
hold on
for i=all_index
    plot(diag(BSRm(:,data_col_real_O(:,i))),rkm_O)
end
title({'BSR sondes';datestr(date_ts_N_O(p_point(1)))})
ylabel('Range (km)')





figure(9985)
bg_o2offm_O = bg_o2off_O.*mask_nodata_O(1,:);
bg_o2offm_O(bg_o2offm_O <= 0) = NaN;
[bg_o2offm_O,TF]=rmmissing(bg_o2offm_O,2);
bg_o2onm_O = bg_o2on_O.*mask_nodata_O(1,:);
bg_o2onm_O(bg_o2onm_O <= 0) = NaN;
bg_o2onm_O=rmmissing(bg_o2onm_O,2);
bg_o2off_molm_O = bg_o2off_mol_O.*mask_nodata_O(1,:);
bg_o2off_molm_O(bg_o2off_molm_O <= 0) = NaN;
bg_o2off_molm_O=rmmissing(bg_o2off_molm_O,2);
bg_o2on_molm_O = bg_o2on_mol_O.*mask_nodata_O(1,:);
bg_o2on_molm_O(bg_o2on_molm_O <= 0) = NaN;
bg_o2on_molm_O=rmmissing(bg_o2on_molm_O,2);

semilogy(datenum(date_ts_N_O(~TF)),bg_o2offm_O)
hold on
semilogy(datenum(date_ts_N_O(~TF)),bg_o2onm_O)
semilogy(datenum(date_ts_N_O(~TF)),bg_o2off_molm_O)
semilogy(datenum(date_ts_N_O(~TF)),bg_o2on_molm_O)
legend('Offline','Online','Offline molecular','Online molecular')
xlim([datenum(date_ts_N_O(1)) datenum(date_ts_N_O(end))]) 
xticks(datenum(date_ts_O(1)): tickSpacing :datenum(date_ts_O(end)))
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
imAlpha=ones(size(T_finalm_O(1:ind_km_max,:)));
imAlpha(isnan(T_finalm_O(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm_O(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(T_finalm_O(1:ind_km_max,:)-1,[],'all'), max(T_finalm_O(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts_N_O),rkm_O(1:ind_km_max),T_finalm_O(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

% for i=1:length(data_col_real_O)
%     if data_col_real_O(i) <= length(date_ts_N_O)
%         xline(datenum(date_ts_N_O(data_col_real_O(i))));
%     end
% end

%Plotting sode locations
for i = 1:size(data_col_real_O,2)
    if ~isnan(data_col_real_O(1,i))
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N_O([data_col_real_O(1,i) data_col_real_O(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
%title(sprintf('Temperature %s to %s (UTC)',date_ts_N_O(1),date_ts_N_O(end)))
%title(sprintf('Temperature retrieval\n%s(UTC)',span_days(1)))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'Temperature (K)')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts_N_O(1)) datenum(date_ts_N_O(end))]) 
xticks(datenum(date_ts_N_O(1)): tickSpacing :datenum(date_ts_N_O(end)))
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
plot(diag(alpha_0_O(:,p_point)),rkm,'-.','LineWidth',2)
hold on
%plot(diag(absorption_O(:,p_point)),rkm)
%%plot(diag(alpha_total_raw_O(:,p_point)),rkm,'LineWidth',2)
%%%plot(absorption_sonde_O{t_index},rkm,'k','LineWidth',2)
plot(diag(alpha_totalm_O(:,p_point)),rkm,'r','LineWidth',2)

%legend('Measured zeroth order absorption','Theoretical absorption from Tg','Pertabative absorption','Alpha 0+1+2','Sonde')%,'theoretical absorption online-offine')
%legend('Measured zeroth order absorption','Alpha 0+1+2','Pertabative absorption','Sonde')%,'theoretical absorption online-offine')
legend('Measured zeroth order absorption','Corrected absorption')
hold off
grid on
ylim([-inf 5])
xlim([0 inf])
title({'O2 absorption profile';datestr(date_ts_N_O(p_point(1)))})
ylabel('Range (km)')
xlabel('Absorption m^{-1}')





figure(8)
plot(rm,o2on_O(:,p_point).*rm_O.^2)
hold on
plot(rm,o2off_O(:,p_point).*rm_O.^2)
%title('Averaged photon returns multiplied by r^2')
title({'Averaged counts returns multiplied by r^2';datestr(date_ts_N_O(p_point))})
legend('Online','Offline')
xlabel('Range m')
ylabel('Counts * r^2')
grid on
hold off

figure(7)
i_range = length(rm_O);
plot(rm,o2on_O(:,p_point))
hold on
plot(rm,o2off_O(:,p_point))
plot(rm,o2on_mol_O(:,p_point))
plot(rm,o2off_mol_O(:,p_point))

plot(rm,o2on_noise_O(1:i_range,p_point),'--')
plot(rm,o2off_noise_O(1:i_range,p_point),'--')
plot(rm,o2on_noise_mol_O(1:i_range,p_point),'--')
plot(rm,o2off_noise_mol_O(1:i_range,p_point),'--')
title({'Averaged counts';datestr(date_ts_N_O(p_point))})
legend('Online','Offline','Online molecular','Offline molecular')
hold off
grid on
xlabel('Range [m]')
ylabel('Counts')


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





% no_cloud_index = [2 4 5 18 22 26 27 28 29 30 31 37 38 39 40 42 43 44 45 ];%4_17 to 4_22
% cloud_index = [1 3 6 7 8 9 10 11 12 13 14 15 16 17 19 20 21 23 24 25 32 33 34 35 36 41];%4_17 to 4_22

figure(6773)
%temp diff sonde
hold on
% index=1;
% for i = 1:length(data_col_real)-1
%     if i==no_cloud_index(index)
%         plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*r')
%         fullTFinal(:,index) = T_finalm(:,data_col_real(i));
%         T_real_final(:,index) = T_real(:,i);
%         index=index+1;
%     else
%         plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*b')
%         fullTFinal(:,index) = T_finalm(:,data_col_real(i));
%         T_real_final(:,index) = T_real(:,i);
%         %index=index+1;
%     end
% 
% end
index=1;
for i = no_cloud_index
    if ~isnan(data_col_real_O(i))
        plot(T_finalm_O(:,data_col_real_O(i)),T_real_O(:,i),'*r')
        fullTFinal_O(:,index) = T_finalm_O(:,data_col_real_O(i));
        T_real_final_O(:,index) = T_real_O(:,i);
        index=index+1;
    end

end
index=1;
for i = cloud_index
    if ~isnan(data_col_real_O(i))
        plot(T_finalm_O(:,data_col_real_O(i)),T_real_O(:,i),'*b')
        fullTFinal_O(:,index) = T_finalm_O(:,data_col_real_O(i));
        T_real_final_O(:,index) = T_real_O(:,i);
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



% figure(6774)
% index =cloud_index;
% for i = 1:length(260:305)
%     for j = 1:length(260:305)
%         
%         tempProb(i,j) = nnz((fullTFinal_O(:,index)<(260+i-1)+0.5)&(fullTFinal_O(:,index)>=(260+i-1)-0.5)& (T_real_final_O(:,index)<(260+j-1)+0.5)&(T_real_final_O(:,index)>=(260+j-1)-0.5));
%     end
% end
% imagesc(260:305,260:305,tempProb')
% colorbar
% set(gca, 'YDir','normal')
% xlabel('Temp retrieval (K)')
% ylabel('Sonde (K)')



% figure(6775)
% 
% %temp diff sonde
% hold on
% index=1;
% for i = cloud_index
%     if ~isnan(data_col_real_O(i))
%         plot(T_finalm_O(:,data_col_real_O(i)),T_real_O(:,i),'*')
%         fullTFinal_O(:,index) = T_finalm_O(:,data_col_real_O(i));
%         T_real_final_O(:,index) = T_real_O(:,i);
%         index=index+1;
%     end
% 
% end
% plot(270:300,270:300)
% plot(273:303,270:300)
% plot(267:297,270:300)
% xlabel('Temperature Retrieval (K)')
% ylabel('Sonde (K)')
% title('Sonde with clouds')
% hold off
% grid on


% figure(6776)
% 
% 
% %temp diff sonde
% hold on
% index=1;
% for i = no_cloud_index
%     if ~isnan(data_col_real_O(i))
%         plot(T_finalm_O(:,data_col_real_O(i)),T_real_O(:,i),'*')
%         fullTFinal_O(:,index) = T_finalm_O(:,data_col_real_O(i));
%         T_real_final_O(:,index) = T_real_O(:,i);
%         index=index+1;
%     end
% 
% end
% plot(270:300,270:300)
% plot(272:302,270:300)
% plot(268:298,270:300)
% xlabel('Temperature Retrieval (K)')
% ylabel('Sonde (K)')
% title('Sonde without clouds')
% hold off
% grid on

figure(7837)
plot(datenum(date_ts_N_O),Ts_O,'Linewidth',2)

xlim([datenum(date_ts_N_O(1)) datenum(date_ts_N_O(end))]) 
xticks(datenum(date_ts_N_O(1)): tickSpacing :datenum(date_ts_N_O(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
ylabel('S Temp (K)')

figure(232)
imagesc(rm,ts,BSR_O,[0 10])
%clim([0 20])
colorbar

%%
alpha_0 = alpha_0_O;
alpha_total_masked = alpha_totalm_O;
alpha_total_raw = alpha_total_raw_O;
date_ts = date_ts_O;
rm = rm_O;
T_final_masked = T_finalm_O;
T_final_raw = T_final_test_O;
ts = ts_O;

save([file 'OwenBoulderData.mat'],'alpha_0','alpha_total_masked','alpha_total_raw','date_ts','rm','T_final_masked','T_final_raw','ts');
