clear all

%9/30/21
%modeling and re-creating Hayman 2020 poisson thinning in matlab
%Owen Cruikshank
path = 'C:\Users\oencr\OneDrive - Montana State University\Research\Papers\Lidar signal procesesing\PoissonLinearProcessing-master\PoissonLinearProcessing-master\data\';
range =double(ncread([path 'mpd05.20181022T12300019921_20181022T15163019921.nc'],'range'));
time=double(ncread([path 'mpd05.20181022T12300019921_20181022T15163019921.nc'],'time'));
Molecular_Counts=ncread([path 'mpd05.20181022T12300019921_20181022T15163019921.nc'],'Molecular_Counts');
Combined_Counts=ncread([path 'mpd05.20181022T12300019921_20181022T15163019921.nc'],'Combined_Counts');



Counts.bg_o2off = mean(Combined_Counts(end-20:end,:));% Take mean of last data points
%Combined_Counts = Combined_Counts - Counts.bg_o2off;       % Background subtracted
%Combined_Counts(Combined_Counts < 0) = 0;         % Minimum of zero

Counts.bg_o2off_mol = mean(Molecular_Counts(end-20:end,:));% Take mean of last data points
%Molecular_Counts = Molecular_Counts - Counts.bg_o2off_mol;       % Background subtracted
%Combined_Counts(Combined_Counts < 0) = 0;         % Minimum of zero

cut = 218;
range = range(4:218);
Molecular_Counts = Molecular_Counts(4:cut,:);
Combined_Counts = Combined_Counts(4:cut,:);

Time.t_step = time(2)-time(1);
Range.rangeBin = range(2)-range(1);
Range.rkm = range/1000;
Time.thr = time/60/60;

Counts.NBins=1;

Counts.o2off_mol = Molecular_Counts;
Counts.o2off = Combined_Counts;

Counts.o2on_mol = Molecular_Counts;
Counts.o2on = Combined_Counts;

[Counts,Ezon,Eton,rangeWidthon,timeWidthon] = poissonThin(Counts);
%%
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
title('off mol')

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

figure(1)
imagesc(Time.thr,Range.rkm,Counts.o2off)
set(gca, 'YDir','normal')
caxis([1 1000])
colorbar
colormap('hot')
set(gca,'ColorScale','log')
title('off comb')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

figure(2)
imagesc(Time.thr,Range.rkm,Counts.o2off_mol)
set(gca, 'YDir','normal')
caxis([1 1000])
colorbar
colormap('hot')
set(gca,'ColorScale','log')
title('off mol')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

figure(11)
imagesc(Time.thr,Range.rkm,Combined_Counts-Counts.bg_o2off)
set(gca, 'YDir','normal')
caxis([1 1000])
colorbar
colormap('hot')
set(gca,'ColorScale','log')
title('off comb')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

figure(22)
imagesc(Time.thr,Range.rkm,Molecular_Counts-Counts.bg_o2off_mol)
set(gca, 'YDir','normal')
caxis([1 1000])
colorbar
colormap('hot')
set(gca,'ColorScale','log')
title('off mol')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

figure(3)
p=500;
plot(Range.rkm,Counts.o2off(:,p)+Counts.bg_o2off(p))
hold on
plot(Range.rkm,Combined_Counts(:,p))
plot(Range.rkm,Counts.o2off_mol(:,p)+Counts.bg_o2off_mol(p))
plot(Range.rkm,Molecular_Counts(:,p))
hold off

figure(4)
semilogx(Counts.Poissonthin.timeWidthoff*Time.t_step,Range.rkm,'.')
hold on
semilogx(Counts.Poissonthin.timeWidthoff_mol*Time.t_step,Range.rkm,'.')
hold off
legend('Comb','Mol')
xlim([1 10^4])
grid on
title('on')
ylabel('Range')

figure(5)
semilogy(Time.thr,Counts.Poissonthin.rangeWidthoff*Range.rangeBin,'.')
hold on
semilogy(Time.thr,Counts.Poissonthin.rangeWidthoff_mol*Range.rangeBin,'.')
hold off
legend('Comb','Mol')
ylim([1 500])
xlim([Time.thr(1) Time.thr(end)])
grid on

