[~,p_point] = min(abs(Options.TimeGrid-10));

overlap_raw = o2on_bgsub./o2on_bgsub_mol;

overlapMinIndex = 20;
%overlap = (o2on./max(o2on(8:end,:),[],1))./(o2on_mol./max(o2on_mol(8:end,:),[],1));
overlap = o2on./o2on_mol;
overlap = overlap./max(overlap(overlapMinIndex:round(end/2),:),[],1);

% Online
overlap_noise = padarray(overlap,[oversample/2,t_avg/2],'replicate');
overlap_filt = filter2(k,overlap_noise,'valid');
overlap = interp2(ts-t_step/2,rm-rangeBin/2,overlap_filt(1:end-1,1:end-1),ts,rm);
overlap = fillmissing(overlap,'nearest',1); % Fill in NaNs in dimension 1
overlap = fillmissing(overlap,'nearest',2); % Fill in NaNs in dimension 2

%o2on_mol = o2on_mol.*overlap(:,p_point).*max(o2on_mol,1)./max(o2on_mol.*overlap(:,p_point),1);

o2on_mol_overlap = o2on_mol .* overlap(:,p_point);

%filter overlap
k = ones(oversample,t_avg)./(oversample*t_avg);     % Kernel

%%
% overlapFile = overlap(:,p_point);
% rmOverlap = rm;
% save('overlapFile4.mat','overlapFile','rmOverlap');

%%
close all
[~,p_point] = min(abs(Options.TimeGrid-9));
[~,p_point2] = min(abs(Options.TimeGrid-10));
[~,p_point3] = min(abs(Options.TimeGrid-11));

figure(15)
subplot(2,1,1)
imagesc(Options.TimeGrid,1:560,DataStructure2.MCS.Channel0.Data')
set(gca, 'YDir','normal')
xline(Options.TimeGrid(p_point),'r');
subplot(2,1,2)
imagesc(Options.TimeGrid,1:560,DataStructure2.MCS.Channel8.Data')
xline(Options.TimeGrid(p_point),'r');


figure(16)
subplot(2,1,1)
plot(DataStructure2.MCS.Channel0.Data(p_point,:))
hold on
plot(DataStructure2.MCS.Channel8.Data(p_point,:))
plot(DataStructure2.MCS.Channel2.Data(p_point,:))
plot(DataStructure2.MCS.Channel10.Data(p_point,:))
legend('0','8','2','10')
ylim([0 2*10^4])
hold off
subplot(2,1,2)
plot(DataStructure2.MCS.Channel0.Data(p_point3,:))
hold on
plot(DataStructure2.MCS.Channel8.Data(p_point3,:))
plot(DataStructure2.MCS.Channel2.Data(p_point3,:))
plot(DataStructure2.MCS.Channel10.Data(p_point3,:))
legend('0','8','2','10')
ylim([0 2*10^4])
hold off

figure(116)
subplot(2,1,1)
plot(rm_raw_o2/1000,o2on_bgsub(:,p_point))
hold on
plot(rm_raw_o2/1000,o2on_bgsub_mol(:,p_point))
plot(rm_raw_o2/1000,o2off_bgsub(:,p_point))
plot(rm_raw_o2/1000,o2off_bgsub_mol(:,p_point))
legend('On com','on mol')
ylim([0 10^4])
subplot(2,1,2)
plot(rm_raw_o2/1000,o2on_bgsub(:,p_point3))
hold on
plot(rm_raw_o2/1000,o2on_bgsub_mol(:,p_point3))
plot(rm_raw_o2/1000,o2off_bgsub(:,p_point3))
plot(rm_raw_o2/1000,o2off_bgsub_mol(:,p_point3))
legend('On com','on mol')
ylim([0 10^4])
hold off

figure(17)
subplot(2,1,1)
plot(rm_raw_o2/1000,o2on_bgsub(:,p_point))
hold on
plot(rm_raw_o2/1000,o2on_bgsub_mol(:,p_point))
plot(rm/1000,o2on(:,p_point))
plot(rm/1000,o2on_mol(:,p_point))
legend('On com','on mol')
ylim([0 10^4])
subplot(2,1,2)
plot(rm_raw_o2/1000,o2on_bgsub(:,p_point3))
hold on
plot(rm_raw_o2/1000,o2on_bgsub_mol(:,p_point3))
plot(rm/1000,o2on(:,p_point3))
plot(rm/1000,o2on_mol(:,p_point3))
legend('On com','on mol')
ylim([0 10^4])
hold off

figure(18)
plot(rm_raw_o2/1000,overlap_raw(:,p_point))
ylim([0 2])

figure(19)
plot(rm/1000,overlap(:,p_point))
hold on
plot(rm/1000,overlap(:,p_point2))
plot(rm/1000,overlap(:,p_point3))
hold off
ylim([0 1.1])
yline(1)
legend('before','after','Location','southeast')
title('Norm com online / norm mol online')

figure(20)
hold on
plot(rm/1000,o2on(:,p_point))
plot(rm/1000,o2on_mol(:,p_point))
plot(rm/1000,o2off(:,p_point))
plot(rm/1000,o2on_mol_overlap(:,p_point))
plot(rm/1000,o2on_mol_overlap(:,p_point)+o2on(:,p_point))
legend('On com','on mol','off','on mol corrected','on com + on mol corrected')
ylim([0 2*10^4])
hold off

figure(21)
hold on
plot(rm/1000,o2on(:,p_point))
plot(rm/1000,o2on_mol(:,p_point))
plot(rm/1000,o2off(:,p_point))
plot(rm/1000,o2off_mol(:,p_point))
legend('On combined','On molecular','Off combined','Off molecular')
ylim([0 1.6*10^4])
xlabel('Range (km)')
ylabel('Counts')
title('Photon counts summed over 60s and filtered over 30 min')
grid on
hold off


figure(22)
hold on
plot(thr,o2on(3:5,:)')
%plot(rm/1000,o2on_mol(:,:))
plot(thr,o2off(3:5,:)')
%plot(rm/1000,o2off_mol(:,:))
%legend('On combined','On molecular','Off combined','Off molecular')
%ylim([0 1.6*10^4])
%xlabel('Range (km)')
%ylabel('Counts')
%title('Photon counts summed over 60s and filtered over 30 min')
grid on
hold off