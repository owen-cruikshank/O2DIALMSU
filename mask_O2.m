function [SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2(o2on,o2off,rm,ts,cloud_p_point,SNR_threshold,SD_threshold,oversample,t_avg)
%File: mask_O2.m
%Date: 03/16/2020
%Author: Owen Cruikshank
%Inputs:
%   -o2on:[photons] (range x time) vector of online return photon counts
%   -o2off:[photons] (range x time) vector of offline return photon counts
%   -rm:[m] (range x 1) vector of range dimension
%   -ts:[s] (1 x time) vector of time dimension
%   -cloud_p_point:[hr] point to plot cloud plots. If 0 will not plot.
%   -SNR_threshold: Signal to noise rato threshold for SNR mask
%   -SD_threshold: Standard deviation threshold for cloud mask
%   -oversample: Number of range bin averaging
%   -t_avg: Number of time bin averaging
%
%Outputs:
%   -SNRm:[none] (range x time) calculated SNR mask
%   -cloud_SDm_above:[none] (range x time) calculated cloud mask
%   -cloud_SDm:[none] (range x time) calculated cloud mask with only clouds
%   -o2on_SNR:[none] (range x time) signal to noise ratio for online counts


% ===================
% === SNR O2 mask ===
% ===================
%SNR_threshold = 30;                      % Set SNR threshold for O2 photon counts
%SNR_threshold = 4;

o2on_SNR = sqrt(o2on);                  % Calculate shot noise SNR
o2on_SNRm = ones(size(o2on_SNR));       % Initialize matrix
o2on_SNRm(o2on_SNR < SNR_threshold & rm > 1000) = 0;% Replace low SNR points with zeros

o2off_SNR = sqrt(o2off);
o2off_SNRm = ones(size(o2off_SNR));
o2off_SNRm(o2off_SNR < SNR_threshold & rm > 1000) = 0;

SNRm = o2off_SNRm .* o2on_SNRm;         %Final combination SNR matrix mask with 0 in places of low snr


% ===================
% === Cloud O2 mask ===
% ===================

%Range corrected returns
o2on_range_corrected = o2on.*rm.^2;
o2off_range_corrected = o2off.*rm.^2;

% log_o2on_range_corrected = log(o2on_range_corrected);
% log_o2off_range_corrected = log(o2off_range_corrected);

%SD_threshold = 6*10^9;

%standard deviation
%stdNeiborhood = true(7,29);                                % NxN window arround each point
stdNeiborhood = true(oversample-1,t_avg-1);                                % NxN window arround each point

o2on_SD = stdfilt(o2on_range_corrected,stdNeiborhood);  % Use a square standard deviation filter on the range corrected return to find a standard deviation
cloud_SDm_on = ones(size(o2on));                        % Initialize mask matrix
%cloud_SDm_on(o2on_SD > 1*10^8) = 0;                       % Set points with a standard deviation greater than the threshold to 0 inicating a cloudy area
cloud_SDm_on(o2on_SD > SD_threshold) = 0;                       % Set points with a standard deviation greater than the threshold to 0 inicating a cloudy area

o2off_SD = stdfilt(o2off_range_corrected,stdNeiborhood);
cloud_SDm_off = ones(size(o2on));
cloud_SDm_off(o2off_SD > SD_threshold) = 0;

cloud_SDm = cloud_SDm_off .* cloud_SDm_on;              % Combine on and offline cloud masks


%Set values above the cloud to -1
cloud_SDm_above = cloud_SDm;
diff_cloud_SDm = diff(cloud_SDm);                           % Calculate difference along range axis
diff_SNRm = diff(SNRm);

[cloud_index_row,cloud_index_col] = find(diff_cloud_SDm);   % Find non-zero values in difference
if ~isempty(cloud_index_row)                                 % Check if any clouds exist
    cloud_index = zeros(1,length(ts));                          % Initalize cloud range index matrix

    index_matrix = ones(length(rm),length(ts)) .* (1:length(rm))'; % Create matrix to match range indecies

    % Set range index value to index of cloud
    cloud_index(cloud_index_col(1))=cloud_index_row(1);         
    for i = 2:length(cloud_index_col)
        if cloud_index_col(i)~=cloud_index_col(i-1)
            cloud_index(cloud_index_col(i))=cloud_index_row(i);
        end
    end
    %cloud_SDm_above(index_matrix > cloud_index - oversample & cloud_index ~= 0)=-1; % Set all values above cloud index to -1
    cloud_SDm_above(index_matrix > cloud_index - oversample & cloud_index ~= 0)=-1; % Set all values above cloud index to -1
    cloud_SDm_above = cloud_SDm_above .* cloud_SDm;                    % Create final matrix where clouds are represended by zeros and all data above represented by -1   
end

[SNR_index_row,SNR_index_col] = find(diff_SNRm);   % Find non-zero values in difference
if ~isempty(SNR_index_row)                                 % Check if any clouds exist
    SNR_index = zeros(1,length(ts));                          % Initalize cloud range index matrix

    index_matrix = ones(length(rm),length(ts)) .* (1:length(rm))'; % Create matrix to match range indecies

    % Set range index value to index of cloud
    SNR_index(SNR_index_col(1))=SNR_index_row(1);         
    for i = 2:length(SNR_index_col)
        if SNR_index_col(i)~=SNR_index_col(i-1)
            SNR_index(SNR_index_col(i))=SNR_index_row(i);
        end
    end
    %cloud_SDm_above(index_matrix > cloud_index - oversample & cloud_index ~= 0)=-1; % Set all values above cloud index to -1
    SNRm(index_matrix > SNR_index )=0; % Set all values above cloud index to -1
                      % Create final matrix where clouds are represended by zeros and all data above represented by -1   
end
    

%---Plotting
if cloud_p_point ~= 0

[~,p_point] = min(abs(cloud_p_point-ts/60/60));%find closest value to 338min for comparison to other program
thr = ts/60/60;
figure(465)
subplot(2,1,1)
imagesc(ts/60/60,rm/1000,cloud_SDm)
hold on
xline(ts(p_point)/60/60,'r');
hold off
set(gca, 'YDir','normal')
pbaspect([3 1 1])
title('Standard deviation mask')
xlabel('Time (UTC hrs)')
ylabel('Range (km)')
subplot(2,1,2)
imagesc(ts/60/60,rm,cloud_SDm_off)
hold on
xline(ts(p_point)/60/60,'r');
hold off
set(gca, 'YDir','normal')
title('o2off standard deviation') 

figure(4465)
imagesc(ts/60/60,rm/1000,cloud_SDm)
set(gca, 'YDir','normal')
%pbaspect([3 1 1])
title('Standard deviation mask')
xlabel('Time (UTC hrs)')
ylabel('Range (km)')


 
figure(565)
plot(rm,o2off_SD(:,p_point))
hold on
plot(rm,o2on_SD(:,p_point))
hold off
legend('on','off')
title('Standard Deviation filter')
xlabel('Range(m)')

figure(665)
plot(rm,o2on_SNR(:,p_point))
hold on
plot(rm,o2off_SNR(:,p_point))
hold off
legend('on','off')
title('Signal to noise ratio')
xlabel('Range(m)')

o2on_SD = o2on_SD.* cloud_SDm_above;
o2on_SD(cloud_SDm_above ~= 1) = NaN;

o2off_SD = o2off_SD.* cloud_SDm_above;
o2off_SD(cloud_SDm_above ~= 1) = NaN;

o2on_masked = o2on_range_corrected.* cloud_SDm_above;
o2on_masked(cloud_SDm_above ~= 1) = NaN;


o2off_masked = o2off_range_corrected.* cloud_SDm_above;
o2off_masked(cloud_SDm_above ~= 1) = NaN;


figure(56)
subplot(2,1,1)
imagesc(ts/60/60,rm,o2on_SD)
hold on
xline(ts(p_point)/60/60,'r');
hold off
set(gca, 'YDir','normal')
title('o2on standard deviation')
subplot(2,1,2)
imagesc(ts/60/60,rm,o2off_SD)
hold on
xline(ts(p_point)/60/60,'r');
hold off
set(gca, 'YDir','normal')
title('o2off standard deviation')


figure(58)
subplot(2,1,1)
imagesc(thr,rm/1000,o2on_range_corrected)
hold on
xline(thr(p_point),'r');
hold off
%caxis([0 1e9])
colorbar
set(gca, 'YDir','normal')
title('Online Counts * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,1,2)
imagesc(thr,rm,o2off_range_corrected)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('o2off range corrected')

figure(558)
imagesc(thr,rm/1000,o2on_range_corrected)
%caxis([0 1e9])
colorbar
%colormap('turbo')
set(gca, 'YDir','normal')
title('Online Counts * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')


figure(588)
subplot(2,1,1)
imagesc(thr,rm/1000,o2on_masked)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')
title('Masked Online Counts * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,1,2)
imagesc(thr,rm,o2off_masked)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
%set(gca,'ColorScale','log')
colorbar
title('o2off range corrected')

figure(5588)

imAlpha=ones(size(o2on_masked));
imAlpha(isnan(o2on_masked))=0;%Set AlphaData to not plot NaNs
%imAlpha(cloud_SDm_above ~= 1)=1;

imagesc(thr,rm/1000,o2on_masked,'AlphaData',imAlpha)
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')


title('Masked Online Counts * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')


figure(59)
subplot(2,1,1)
imagesc(thr,rm,o2on)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('o2on')
subplot(2,1,2)
imagesc(thr,rm,o2off)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('o2off')

figure(60)
subplot(2,1,1)
imagesc(thr,rm,o2on_SNRm)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('o2on_SNRm')
subplot(2,1,2)
imagesc(thr,rm,o2off_SNRm)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('oo2off_SNRm')
end

end