load('330T.mat')
T330 = Temperature;
load('240T.mat')
T240 = Temperature;

load('330T2.mat')
T330 = Temperature;
load('240T2.mat')
T240 = Temperature;
%load('normalT.mat')

figure()
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
r_min_plot = .5; % Max range to plot [km]
[~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
%imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
%imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
%imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
%colorLimits = [min(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1
colorLimits =[-10 10];
imagesc(datenum(Time.date_ts),Range.rkm(1:end),T240.T_final_test(1:end,:)-T330.T_final_test(1:end,:),colorLimits) %Plot
hold on
%sonde lines
 %plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% % %yline(.6)
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% % end
% % xline(738397.94)
% % xline(738398.17)
% % xline(738400.65)
% % xline(738401.18)
colors = colormap; % Set colors to colormap
%colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colors = [1 0 0
    1 1 1
    0 0 1];

colors = [[linspace(0,1,20) linspace(1,1,20)]' [linspace(0,1,20) linspace(1,0,20)]' [linspace(1,1,20) linspace(1,0,20)]'];

colormap(colors); % Make colormap new colors
set(gca, 'YDir','normal') % Set y axis increasing
%set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
xticks(datenum(Time.date_ts(1)):Format.tickSpacing:datenum(Time.date_ts(end)))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')
title(sprintf('Temperature retrieval')) %Add title
xlabel('Time (UTC hours)')
ylabel('Range (km)')
title(colorbar,'Temperature (K)') % Add title to colorbar
hold off