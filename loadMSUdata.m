function [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data,Options] = loadMSUdata(span_days,Options,Constant)
%====================
%==== Reading Files =====
%====================
disp('Reading in files')
maindirectory =pwd;
cd ../
Options.path = fullfile(pwd,'Data','MSU data','RSync','NetCDFOutput'); %Path for instument data
Options.weatherPath = fullfile(pwd, 'Data'); %path for weather station data
Options.sondepath = fullfile(pwd ,'Data','MSU data','Radiosondes'); %path for radiosonde data

cd( maindirectory) %back to main directory

Options.MPDname = 'MSU';
Options.BinTotal = 560;
%Options.BinTotal = 400;
%Options.BinTotal = 490;
%Load raw data from NetCDF files
[Data, Options] = loadMSUNETcdf(span_days,Options);


disp('Creating time and range vectors')
% ====================
% === Time Vectors ===
% ====================
Time.ts = Data.MCS.Channel0.NewTimeGrid*60*60; %(s) Time vector from data
Time.t_step = Time.ts(2)-Time.ts(1); %(s) time bin size
Time.i_time = length(Time.ts); % length of time vector
Time.date_ts = span_days(1) + seconds(Time.ts); %(datetime) time vector in datetime format
Time.thr = Time.ts/60/60; %(hr)


% =====================
% === Range Vectors ===
% =====================
% Create range vector from length of data bins
Range.nsPerBin = 250; %[ns] bin length in nanosections
Range.NBins = floor(Options.BinTotal/Options.intRange); %number of range bins in vector
Range.rangeBin = (Constant.c * Range.nsPerBin(1)*10^-9)/2; %(m)range bin length

Range.rm_raw_o2 = -150:Range.rangeBin:Range.NBins(1)*Range.rangeBin-150-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = -150-30:Range.rangeBin:Range.NBins(1)*Range.rangeBin-150-30-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = -150-30/2+Range.rangeBin:Range.rangeBin:Range.NBins(1)*Range.rangeBin-150-30/2-Range.rangeBin+Range.rangeBin;    %[m] Create range vector
%%%%%%%%%%%%%%%%%%%%%%%
Range.rm_raw_o2 = -75-30/2:Range.rangeBin:Range.NBins(1)*Range.rangeBin+-75-30/2-Range.rangeBin;    %[m] Create range vector
%%%%%%%%%%%%%%%%%%%%%%%%
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector
Range.r_max = 6000;                                       %[m] Max range 
Range.rm = Range.rm_raw_o2(Range.rm_raw_o2<=Range.r_max & Range.rm_raw_o2>0);     %[m] Shorten range vector to max range

%=== Integrate range vector ==
Range.rm = Range.rm(1:2:end);%integrate to new
Range.rangeBin = Range.rangeBin*2;

Range.i_range = length(Range.rm);                               %[none] Size of range vector
Range.rkm = Range.rm/1000;


%==============================
%== MSU Weather station data ==
%==============================
disp('Loading weather Station Data')
[weather_Temperature_interp, weather_absPressure_interp, ~] = ORSLweatherv3(span_days,Time.ts,Options.weatherPath);
%Use for weather from weather underground data
%[weather_Temperature_interp, weather_absPressure_interp, weather_WV_interp] = wunderWeather(span_days,ts,path);
%Use if there is no data from weather station
%weather_Temperature_interp = 272.310000000000*ones(1,length(ts)) -273.15;
%weather_absPressure_interp = 0.836910930175179*ones(1,length(ts)).* 1013.25;
%weather_WV_interp = zeros(1,Time.i_time);


disp('Calculating model')
% === Calculating temperature and pressure model ===
Model.Ts = weather_Temperature_interp + 273.15 ;          %surface temperature from weather station [K]
Model.Ps = weather_absPressure_interp / 1013.25;         %absolute surface pressure from weather station [atm]
lapseRate = -6.5;                                   %[K/km] Guess adiabatic lapse rate  typically -6.5 up to 10km
lapseRate = -10;    
lapseRate = lapseRate / 1000;                       %[K/m] 

Model.T = Model.Ts + lapseRate .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
Model.P = Model.Ps .* (Model.Ts./Model.T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  

Model.lapseRate = lapseRate;
% == water vapor model ==
%Pws = exp(77.3450+0.0057.*Model.T-7235./Model.T)./Model.T.^8.2;
%WV = weather_WV_interp.*0.0022.*Pws./T./100;
%Model.WV = weather_WV_interp.*Pws./Model.T./100/Constant.kb;%water vapror molc/m^3
Model.WV = zeros(size(Model.T));
%%%%%%%WV = ones(i_range,i_time).*WV_sonde(:,1);

%=============================
%== Cobleigh sonde =========
%===========================
disp('Loading Sonde data')
[sonde_datetime,sondeStruc] =  COBradiosonde(Options.sondepath,span_days);
for i = 1:numel(sonde_datetime) % Loop over number of sondes in time period
    if isdatetime(sonde_datetime(i)) %== Check if sonde exists
        % ===Subtract first range value (site elevation) from whole vector
        %rm_sgp{i} = sondeStruc(i).Height - sondeStruc(i).Height(1);
        rm_sgp{i} = sondeStruc(i).Height - 1524;
        %===convert to same units====
        sondeStruc(i).P = sondeStruc(i).P./1013.25;%atm
        % ==Collect radiosonde surface measurements==
        T_sgp_surf(i) = sondeStruc(i).T(1);
        P_sgp_surf(i) = sondeStruc(i).P(1);
        % ==Custom interpolation function==
        [T_sonde_int{i},P_sonde_int{i},WV_sonde_int{i},rm_sonde_int{i}] = interp_sonde2(sondeStruc(i).T,sondeStruc(i).P,sondeStruc(i).WV,rm_sgp{i},Range.rangeBin);  
        
        [rm_sgp{i},IA,~] = unique(rm_sgp{i});
        sondeStruc(i).T = sondeStruc(i).T(IA);
        sondeStruc(i).P = sondeStruc(i).P(IA);
        sondeStruc(i).WV = sondeStruc(i).WV(IA);
        T_sonde_int{i} = interp1(rm_sgp{i},sondeStruc(i).T,Range.rm,'makima',nan);
        P_sonde_int{i} = interp1(rm_sgp{i},sondeStruc(i).P,Range.rm,'makima',nan);
        WV_sonde_int{i} = interp1(rm_sgp{i},sondeStruc(i).WV,Range.rm,'makima',nan);
        
        if length(T_sonde_int{i})<Range.i_range % ==If sonde does not reach full lidar range
            disp('ran')
            Sonde.T_sonde(1:length(T_sonde_int{i}),i) = T_sonde_int{i};
            Sonde.T_sonde(length(T_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(T_sonde_int{i}),1);
            Sonde.P_sonde(1:length(T_sonde_int{i}),i) = P_sonde_int{i};
            Sonde.P_sonde(length(P_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(P_sonde_int{i}),1);
            Sonde.WV_sonde(1:length(T_sonde_int{i}),i) = WV_sonde_int{i};
            Sonde.WV_sonde(length(WV_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(WV_sonde_int{i}),1);
            Sonde.Tsurf(:,i) = sondeStruc(i).T(1);
            Sonde.Psurf(:,i) = sondeStruc(i).P(1);
            Sonde.AbsHum(:,i) = Sonde.WV_sonde(:,i).*Constant.mWV*1000; %[g/m3]
        else
            Sonde.T_sonde(:,i) = T_sonde_int{i}(1:Range.i_range);
            Sonde.P_sonde(:,i) = P_sonde_int{i}(1:Range.i_range);
            Sonde.WV_sonde(:,i) = WV_sonde_int{i}(1:Range.i_range);
            Sonde.Tsurf(:,i) = sondeStruc(i).T(1);
            Sonde.Psurf(:,i) = sondeStruc(i).P(1);
            Sonde.AbsHum(:,i) = Sonde.WV_sonde(:,i).*Constant.mWV*1000; %[g/m3]
        end
        %===interp sonde time
%         %[rm_sgp{i},IA,~] = unique(rm_sgp{i});
%         sonde_time(1:length(rm_sonde_int{i}),i) = interp1(rm_sgp{i},sondeStruc(i).time(IA),rm_sonde_int{i})';  
%         if length(sonde_time) < Range.i_range
%             sonde_time = [sonde_time; sonde_time(end).*ones(Range.i_range-length(sonde_time),1)];
%         end
%          %Find index of sonde in time vector
%          for j = 1:Range.i_range
%             [~, Sonde.sonde_ind(j,i)]=min(abs(sonde_datetime(i)+seconds(sonde_time(j,i))-Time.date_ts));
%          end

        sonde_time(1:length(rm_sonde_int{i})) = interp1(rm_sgp{i},sondeStruc(i).time(IA),rm_sonde_int{i})';  
        if length(sonde_time) < Range.i_range
            sonde_time = [sonde_time; sonde_time(end).*ones(Range.i_range-length(sonde_time))];
        end
         %Find index of sonde in time vector
         for j = 1:Range.i_range
            [~, Sonde.sonde_ind(j,i)]=min(abs(sonde_datetime(i)+seconds(sonde_time(j))-Time.date_ts));
         end

    else
        Sonde.sonde_ind = [];
    end 
end

%=== set model WV to sonde interpolation
if ~isempty(Sonde.sonde_ind)
    if length(Sonde.sonde_ind(1,:))==1 %==if only one sonde
        Model.WV = ones(size(Model.WV)).*Sonde.WV_sonde;
        Model.WV = fillmissing(Model.WV,'linear',1);
        Model.WV = fillmissing(Model.WV,'linear',2);
    else
        Model.WV = interp2(Time.ts(Sonde.sonde_ind(1,:)),Range.rm,Sonde.WV_sonde,Time.ts,Range.rm,'linear');
        Model.WV = fillmissing(Model.WV,'linear',1);
        Model.WV = fillmissing(Model.WV,'linear',2);
    end
end

%%%%Model.WV = zeros(size(Model.WV));
%%
% %==== After pulse correction ====
% load('AfterPulse2.mat','PulseOn','PulseOff','PulseOnMol','PulseOffMol')
% Data.MCS.Channel2.Data = Data.MCS.Channel2.Data - PulseOn+1;
% Data.MCS.Channel10.Data = Data.MCS.Channel10.Data - PulseOff+1;
% Data.MCS.Channel0.Data = Data.MCS.Channel0.Data - PulseOnMol+1;
% Data.MCS.Channel8.Data = Data.MCS.Channel8.Data - PulseOffMol+1;

% 
% load('C:\Users\Owen C\Downloads\Afterpulsing_correction_08082022.mat')
% %load('C:\Users\Owen\Downloads\Afterpulsing_correction_08082022.mat')
% % Correction_Nc_off =[Correction_Nc_off'; zeros(560-490,1)];
% % Correction_Nc_on =[Correction_Nc_on'; zeros(560-490,1)];
% % Correction_Nm_on = [Correction_Nm_on'; zeros(560-490,1)];
% % Correction_Nm_off = [Correction_Nm_off'; zeros(560-490,1)];
% 
% Correction_Nc_off =Correction_Nc_off';
% Correction_Nc_on =Correction_Nc_on'; 
% Correction_Nm_on = Correction_Nm_on';
% Correction_Nm_off = Correction_Nm_off'; 
% 
% Data.MCS.Channel10.Data = Data.MCS.Channel10.Data-Correction_Nc_off;
% Data.MCS.Channel2.Data = Data.MCS.Channel2.Data-Correction_Nc_on;
% Data.MCS.Channel0.Data = Data.MCS.Channel0.Data-Correction_Nm_on;
% Data.MCS.Channel8.Data = Data.MCS.Channel8.Data-Correction_Nm_off;
%%
%=========================
%=== Calculate background
%=========================
% % Counts.bg_o2off = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2off_bgsub = Data.MCS.Channel10.Data - Counts.bg_o2off;       % Background subtracted
% % Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
% % 
% % Counts.bg_o2on = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2on_bgsub = Data.MCS.Channel2.Data - Counts.bg_o2on;       % Background subtracted
% % Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
% % 
% % Counts.bg_o2on_mol = mean(Data.MCS.Channel0.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2on_bgsub_mol = Data.MCS.Channel0.Data - Counts.bg_o2on_mol;       % Background subtracted
% % Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
% % 
% % Counts.bg_o2off_mol = mean(Data.MCS.Channel8.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2off_bgsub_mol = Data.MCS.Channel8.Data - Counts.bg_o2off_mol;       % Background subtracted
% % Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero


% Counts.bg_o2off = mean(Data.MCS.Channel10.Data(400-40:400,:));% Take mean of last data points
% Counts.o2off_bgsub = Data.MCS.Channel10.Data - Counts.bg_o2off;       % Background subtracted
% Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2on = mean(Data.MCS.Channel2.Data(400-40:400,:));% Take mean of last data points
% Counts.o2on_bgsub = Data.MCS.Channel2.Data - Counts.bg_o2on;       % Background subtracted
% Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2on_mol = mean(Data.MCS.Channel0.Data(400-40:400,:));% Take mean of last data points
% Counts.o2on_bgsub_mol = Data.MCS.Channel0.Data - Counts.bg_o2on_mol;       % Background subtracted
% Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2off_mol = mean(Data.MCS.Channel8.Data(400-40:400,:));% Take mean of last data points
% Counts.o2off_bgsub_mol = Data.MCS.Channel8.Data - Counts.bg_o2off_mol;       % Background subtracted
% Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero

Counts.bg_o2off = mean(Data.MCS.Channel10.Data(490-40:490,:));% Take mean of last data points
Counts.o2off_bgsub = Data.MCS.Channel10.Data - Counts.bg_o2off;       % Background subtracted
Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero

Counts.bg_o2on = mean(Data.MCS.Channel2.Data(490-40:490,:));% Take mean of last data points
Counts.o2on_bgsub = Data.MCS.Channel2.Data - Counts.bg_o2on;       % Background subtracted
Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero

Counts.bg_o2on_mol = mean(Data.MCS.Channel0.Data(490-40:490,:));% Take mean of last data points
Counts.o2on_bgsub_mol = Data.MCS.Channel0.Data - Counts.bg_o2on_mol;       % Background subtracted
Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero

Counts.bg_o2off_mol = mean(Data.MCS.Channel8.Data(490-40:490,:));% Take mean of last data points
Counts.o2off_bgsub_mol = Data.MCS.Channel8.Data - Counts.bg_o2off_mol;       % Background subtracted
Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero

%%

% ========integrate to new range
inc = 1;
ii = 1:2:length(Range.rm_raw_o2);
o2on_intp2 = zeros(length(ii),length(Counts.o2on_bgsub(1,:)));
o2off_intp2 = zeros(length(ii),length(Counts.o2on_bgsub(1,:)));
o2on_intp2_mol = zeros(length(ii),length(Counts.o2on_bgsub(1,:)));
o2off_intp2_mol = zeros(length(ii),length(Counts.o2on_bgsub(1,:)));
for ii = 1:2:length(Range.rm_raw_o2)
    o2on_intp2(inc,:) = (Counts.o2on_bgsub(ii,:)+Counts.o2on_bgsub(ii+1,:))/2;
    o2off_intp2(inc,:) = (Counts.o2off_bgsub(ii,:)+Counts.o2off_bgsub(ii+1,:))/2;
    o2on_intp2_mol(inc,:) = (Counts.o2on_bgsub_mol(ii,:)+Counts.o2on_bgsub_mol(ii+1,:))/2;
    o2off_intp2_mol(inc,:) = (Counts.o2off_bgsub_mol(ii,:)+Counts.o2off_bgsub_mol(ii+1,:))/2;
    inc = inc+1;
end
Counts.o2on_bgsub = o2on_intp2;
Counts.o2off_bgsub = o2off_intp2;
Counts.o2on_bgsub_mol = o2on_intp2_mol;
Counts.o2off_bgsub_mol = o2off_intp2_mol;

Range.rm_raw_o2 = Range.rm_raw_o2(1:2:end)+Range.rangeBin./2;

%%

% ===Interpolating to shorter range vector===
Counts.o2on_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub,Time.ts,Range.rm);
Counts.o2on_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub_mol,Time.ts,Range.rm);
Counts.o2off_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub,Time.ts,Range.rm);
Counts.o2off_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub_mol,Time.ts,Range.rm);

%integrate Bins to new range
Counts.NBins = Data.MCS.Channel0.NBins*2;
%Counts.NBins = Data.MCS.Channel0.NBins;


%%
%Dead time correction
% % % deadTime = 22e-9; %SPCM-AQRH-13 dead time
% % % Counts.o2onCR = Counts.o2on_noise.*Counts.NBins.*250e-9.*14000;%Count rate, Counts*NuberTimeSummedbins*Length of bin(250ns)*profiles per histogram
% % % Counts.o2on_noise = Counts.o2on_noise ./(1-(deadTime.*Counts.o2onCR));
% % % Counts.o2offCR = Counts.o2off_noise.*Counts.NBins.*250e-9.*14000;%Count rate, Counts*NuberTimeSummedbins*Length of bin(250ns)*profiles per histogram
% % % Counts.o2off_noise = Counts.o2off_noise ./(1-(deadTime.*Counts.o2offCR));
% % % Counts.o2on_molCR = Counts.o2on_noise_mol.*Counts.NBins.*250e-9.*14000;%Count rate, Counts*NuberTimeSummedbins*Length of bin(250ns)*profiles per histogram
% % % Counts.o2on_noise_mol = Counts.o2on_noise_mol ./(1-(deadTime.*Counts.o2on_molCR));
% % % Counts.o2off_molCR = Counts.o2off_noise_mol.*Counts.NBins.*250e-9.*14000;%Count rate, Counts*NuberTimeSummedbins*Length of bin(250ns)*profiles per histogram
% % % Counts.o2off_noise_mol = Counts.o2off_noise_mol ./(1-(deadTime.*Counts.o2off_molCR));
%%
%===== Afterpulse Correction =====
% load('AfterPulse.mat','pulseON','pulseOFF','pulseON_mol','pulseOFF_mol')
% Counts.o2on_noise(5:end,:) = Counts.o2on_noise(5:end,:)-pulseON(5:end,:)*1;
% Counts.o2off_noise(5:end,:) = Counts.o2off_noise(5:end,:)-pulseON(5:end,:)*1;
% Counts.o2on_noise_mol(5:end,:) = Counts.o2on_noise_mol(5:end,:)-pulseON_mol(5:end,:)*1;
% Counts.o2off_noise_mol(5:end,:) = Counts.o2off_noise_mol(5:end,:)-pulseOFF_mol(5:end,:)*1;
% 
% Counts.o2on_noise(Counts.o2on_noise<0)=0;
% Counts.o2off_noise(Counts.o2off_noise<0)=0;
% Counts.o2on_noise_mol(Counts.o2on_noise_mol<0)=0;
% Counts.o2off_noise_mol(Counts.o2off_noise_mol<0)=0;




%%
%=====Create Spectrum vectors=====

%lambda_online = interp1(Options.TimeGrid,Data.Laser.O2Online.WavelengthActual,Time.ts/60/60);
%lambda_offline = interp1(Options.TimeGrid,Data.Laser.O2Offline.WavelengthActual,Time.ts/60/60);

Spectrum.lambda_online = 769.7958 *ones(size(Time.ts));
Spectrum.lambda_offline = 770.1085 *ones(size(Time.ts));

% Spectrum.lambda_wvon = fillmissing(filloutliers(Data.Laser.WVOnline.WavelengthActual,'linear','movmedian',5),'linear');
% Spectrum.lambda_wvoff = fillmissing(filloutliers(Data.Laser.WVOffline.WavelengthActual,'linear','movmedian',5),'linear');

Spectrum.lambda_wvon = 828.1959 *ones(size(Time.ts));
Spectrum.lambda_wvoff = 828.2951 *ones(size(Time.ts));


Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber

Spectrum.nu_wvon = 10^7./Spectrum.lambda_wvon;                    %[cm-1] Online wavenumber
Spectrum.nu_wvoff = 10^7./Spectrum.lambda_wvon;                  %[cm-1] Offline wavenumber

nuMin = Spectrum.nu_online-0.334;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nuwvMin = mean(Spectrum.nu_wvon)-0.334;                                 %[cm-1] Scan lower bound
nuwvMax = mean(Spectrum.nu_wvon)+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
nu_scanwv = (nuwvMin:Spectrum.nuBin:nuwvMax);                      %[cm-1](1 x nu) Scan vector

nuMin_off = Spectrum.nu_offline-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334;                                 %[cm-1] Scan upper bound
nu_scan_off = (nuMin_off:Spectrum.nuBin:nuMax_off);

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension

Spectrum.nu_scanwv_3D_short = permute(nu_scanwv, [3 1 2]);       %[cm-1] putting scan in third dimension

Spectrum.lambda_scan_3D_short = 10^7./Spectrum.nu_scan_3D_short;
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);         %[none] length of scan vector

Spectrum.lambda_scanwv_3D_short = 10^7./Spectrum.nu_scan_3D_short;

Spectrum.del_nu = Spectrum.nu_scan_3D_short-Spectrum.nu_online;                %[1/cm] difference from center
Spectrum.del_lambda = Spectrum.lambda_scan_3D_short-Spectrum.lambda_online;

[~,Spectrum.online_index] = min(abs(Spectrum.nu_online - Spectrum.nu_scan_3D_short),[],3);%finding index of online wavenumber
[~,Spectrum.offline_index] = min(abs(Spectrum.nu_offline - Spectrum.nu_scan_3D_short_off),[],3);%finding index of online wavenumber

[~,Spectrum.online_indexwv] = min(abs(Spectrum.nu_wvon - Spectrum.nu_scanwv_3D_short),[],3);%finding index of online wavenumber


%%
%===== Calculate Model absorption from Model T and P =======
%%%%Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption

Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV);

Model.absorption_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption
Model.transmission = exp(-cumtrapz(Range.rm,Model.absorption));

%===== Calucation Model absorption for radiosondes =======
for i=1:numel(sonde_datetime) 
        if isdatetime(sonde_datetime(i)) %Check if there are any sondes
            Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(Sonde.sonde_ind(:,i)),Model.WV(:,Sonde.sonde_ind(:,i)))); %[m-1] Funcrtion to calculate theoretical absorption
             Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(Sonde.sonde_ind(:,i)),Sonde.WV_sonde(:,i))); %[m-1] Funcrtion to calculate theoretical absorption
            Sonde.trasmission_sonde{i} = exp(-cumtrapz(Range.rm,Sonde.absorption_sonde{i})); %O2 transmission
        else
            Sonde.absorption_sonde{i} = nan(Range.i_range,1);
            Sonde.trasmission_sonde{i} = nan(Range.i_range,1);
            Sonde.T_sonde = nan(Range.i_range,1);
            Sonde.P_sonde = nan(Range.i_range,1);
            Sonde.WV_sonde = nan(Range.i_range,1);
        end
end

%%
%====Count filtering====
k = ones(Options.oversample,Options.t_avg)./(Options.oversample*Options.t_avg);     % Kernel
Counts.o2on = filter2(k,Counts.o2on_noise,'same');
Counts.o2on_mol = filter2(k,Counts.o2on_noise_mol,'same');
Counts.o2off = filter2(k,Counts.o2off_noise,'same');
Counts.o2off_mol = filter2(k,Counts.o2off_noise_mol,'same');

%%
%Afterpulse
% load('AfterPulse.mat')
% Counts.o2on = Counts.o2on + pulseON(1:2:160);
% Counts.o2off = Counts.o2off + pulseON(1:2:160);

%%
Counts.wvon = nan(size(Counts.o2on));
Counts.wvoff = nan(size(Counts.o2on));
Counts.wvon_noise = nan(size(Counts.o2on_noise));
Counts.wvoff_noise = nan(size(Counts.o2on_noise));
Counts.bg_wvon = nan(size(Counts.bg_o2on));
Counts.bg_wvoff = nan(size(Counts.bg_o2on));

%%
% for iii = 1:Time.i_time
% Counts.deconvo2on(:,iii) =deconv(Counts.o2on(:,iii),[0.5 0.5]);
% Counts.deconvo2off(:,iii) =deconv(Counts.o2off(:,iii),[0.5 0.5]);
% Counts.deconvo2on_mol(:,iii) =deconv(Counts.o2on_mol(:,iii),[0.5 0.5]);
% Counts.deconvo2off_mol(:,iii) =deconv(Counts.o2off_mol(:,iii),[0.5 0.5]);
% end
% ending = Counts.deconvo2on(end-1:end,:);
% endingoff = Counts.deconvo2off(end-1:end,:);
% 
% ending_mol = Counts.deconvo2on_mol(end-1:end,:);
% endingoff_mol = Counts.deconvo2off_mol(end-1:end,:);
% for iii = 1:2:Range.i_range-2
%     Counts.deconvo2on(iii:iii+1,:) = Counts.deconvo2on(iii:iii+1,:)+ending;
%     Counts.deconvo2off(iii:iii+1,:) = Counts.deconvo2off(iii:iii+1,:)+endingoff;
%     Counts.deconvo2on_mol(iii:iii+1,:) = Counts.deconvo2on_mol(iii:iii+1,:)+ending_mol;
%     Counts.deconvo2off_mol(iii:iii+1,:) = Counts.deconvo2off_mol(iii:iii+1,:)+endingoff_mol;
% end
% 
% Counts.deconvo2on(end+1:end+2,:) = Counts.deconvo2on(end-1:end,:);
% Counts.deconvo2off(end+1:end+2,:) = Counts.deconvo2off(end-1:end,:);
% Counts.deconvo2on_mol(end+1:end+2,:) = Counts.deconvo2on_mol(end-1:end,:);
% Counts.deconvo2off_mol(end+1:end+2,:) = Counts.deconvo2off_mol(end-1:end,:);

% Counts.o2on = Counts.deconvo2on(1:end-1,:);
% Counts.o2off =Counts.deconvo2off(1:end-1,:);
% Counts.o2on_mol = Counts.deconvo2on_mol(1:end-1,:);
% Counts.o2off_mol = Counts.deconvo2off_mol(1:end-1,:);


%%
%====== Calucate any appy optimal filtering based on Poisson thinning ====
%Counts = poissonThin(Counts);

Counts.Poissonthin.timeWidthon = nan(size(Counts.o2on,1));
Counts.Poissonthin.timeWidthoff = nan(size(Counts.o2on,1));
Counts.Poissonthin.timeWidthon_mol = nan(size(Counts.o2on,1));
Counts.Poissonthin.timeWidthoff_mol = nan(size(Counts.o2on,1));
Counts.Poissonthin.rangeWidthon = nan(size(Counts.o2on,2));
Counts.Poissonthin.rangeWidthoff = nan(size(Counts.o2on,2));
Counts.Poissonthin.rangeWidthon_mol = nan(size(Counts.o2on,2));
Counts.Poissonthin.rangeWidthoff_mol = nan(size(Counts.o2on,2));

Counts.Poissonthin.timeWidthwvon = nan(size(Counts.o2on,1));
Counts.Poissonthin.timeWidthwvoff = nan(size(Counts.o2on,1));
Counts.Poissonthin.rangeWidthwvon = nan(size(Counts.o2on,2));
Counts.Poissonthin.rangeWidthwvoff = nan(size(Counts.o2on,2));

 Counts.foff=nan(size(Counts.o2on));
 Counts.foff_mol=nan(size(Counts.o2on));

%%
% Deconvolution estimate
% pulse = [0.5;0.5];
% onConv = conv2(Counts.o2on,pulse);
% Counts.o2on = 2*Counts.o2on - onConv(2:end,:);
% Counts.o2on(Counts.o2on<0)=0;
% offConv = conv2(Counts.o2off,pulse);
% Counts.o2off = 2*Counts.o2off - offConv(2:end,:);
% Counts.o2off(Counts.o2off<0)=0;
% 
% on_molConv = conv2(Counts.o2on_mol,pulse);
% Counts.o2on_mol = 2*Counts.o2on_mol - on_molConv(2:end,:);
% Counts.o2on_mol(Counts.o2on_mol<0)=0;
% off_molConv = conv2(Counts.o2off_mol,pulse);
% Counts.o2off_mol = 2*Counts.o2off_mol - off_molConv(2:end,:);
% Counts.o2off_mol(Counts.o2off_mol<0)=0;

%%
%====Overlap Correction ======
% % load('overlapCorrMSU6_18_21.mat','overlapCorr')
% load('overlapCorrMSU6_21_21.mat','overlapCorr')
% %Counts.o2off = Counts.o2off.*overlapCorr;

% load('overlapCorrMSU6_24_21.mat','overlapCorr')
% load('overlapCorrMSU6_24_21_1000.mat','overlapCorr')
%  [~,constAlt]=min(abs(2000-Range.rm));
%  Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*overlapCorr(1:constAlt,:);

%%
%=====================
%= Backscatter ratio =
%=====================
if span_days(1)<datetime(2020,10,6,'TimeZone','UTC')
    load('overlap.mat','Correction')
    overlapcorrection = interp1(rm_raw_o2,Correction,rm);
    [HSRL.Bm,HSRL.Ba,HSRL.BSR]= BackscatterRatioV3(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);

    %????????

%elseif span_days(1) >= datetime(2022,5,11,'TimeZone','UTC')
    elseif span_days(1) >= datetime(2022,4,12,'TimeZone','UTC')

    Atmosphere.Pressure = Model.P./0.009869233; 
    Atmosphere.Temperature = Model.T;
    Counts.Nc_on = Counts.o2on;
    Counts.Nc_off = Counts.o2off;
    Counts.Nm_on = Counts.o2on_mol;
    Counts.Nm_off = Counts.o2off_mol;
    Options.t_step = 1;
    [HSRL] = HSRL_retrieval(Counts,Atmosphere,Options);

HSRL.BSRf = nan(size(HSRL.BSR));

        HSRL.Bm828 = HSRL.Bm *(770/828)^4;
    HSRL.Ba828 = HSRL.Ba*(770/828);
    HSRL.BSR828 = HSRL.Ba828./HSRL.Bm828+1;


elseif span_days(1)>=datetime(2022,4,11,'TimeZone','UTC')
    LidarData.Range = Range.rm;
    LidarData.Time = Time.ts;
    LidarData.OfflineCombinedTotalCounts = Counts.o2off;
    LidarData.OfflineMolecularTotalCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    [LidarData]=BackscatterRetrieval05202022(LidarData,WeatherData);
    HSRL.BSR = LidarData.UnmaskedBackscatterRatio;
    HSRL.Ba = LidarData.UnmaskedAerosolBackscatterCoefficient;
    HSRL.Bm = LidarData.MolecularBackscatterCoefficient;

    HSRL.Bm828 = LidarData.MolecularBackscatterCoefficient *(770/828)^4;
    HSRL.Ba828 = LidarData.UnmaskedAerosolBackscatterCoefficient*(770/828);
    HSRL.BSR828 = HSRL.Ba828./HSRL.Bm828+1;

%     LidarData.OfflineCombinedTotalCounts = Counts.foff;
%     LidarData.OfflineMolecularTotalCounts = Counts.foff_mol;
%     [LidarData]=BackscatterRetrieval03112022(LidarData,WeatherData);
%    HSRL.BSRf = LidarData.UnmaskedBackscatterRatio;
    HSRL.BSRf = nan(size(HSRL.BSR));



elseif span_days(1)>=datetime(2022,3,11,'TimeZone','UTC')
    LidarData.Range = Range.rm;
    LidarData.Time = Time.ts;
    LidarData.OfflineCombinedTotalCounts = Counts.o2off;
    LidarData.OfflineMolecularTotalCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    [LidarData]=BackscatterRetrieval03112022(LidarData,WeatherData);
    HSRL.BSR = LidarData.UnmaskedBackscatterRatio;
    HSRL.Ba = LidarData.UnmaskedAerosolBackscatterCoefficient;
    HSRL.Bm = LidarData.MolecularBackscatterCoefficient;

    HSRL.Bm828 = LidarData.MolecularBackscatterCoefficient *(770/828)^4;
    HSRL.Ba828 = LidarData.UnmaskedAerosolBackscatterCoefficient*(770/828);
    HSRL.BSR828 = HSRL.Ba828./HSRL.Bm828+1;

    LidarData.OfflineCombinedTotalCounts = Counts.foff;
    LidarData.OfflineMolecularTotalCounts = Counts.foff_mol;
    [LidarData]=BackscatterRetrieval03112022(LidarData,WeatherData);
    HSRL.BSRf = LidarData.UnmaskedBackscatterRatio;

elseif span_days(1)>=datetime(2021,7,6,'TimeZone','UTC')
    LidarData.Range = Range.rm;
    LidarData.Time = Time.ts;
    LidarData.OfflineCombinedTotalCounts = Counts.o2off;
    LidarData.OfflineMolecularTotalCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    [LidarData]=BackscatterRetrievalRayleighBrillouin070621(LidarData,WeatherData);
    HSRL.BSR = LidarData.UnmaskedBackscatterRatio;
    HSRL.Ba = LidarData.UnmaskedAerosolBackscatterCoefficient;
    HSRL.Bm = LidarData.MolecularBackscatterCoefficient;

    HSRL.Bm828 = LidarData.MolecularBackscatterCoefficient *(770/828)^4;
    HSRL.Ba828 = LidarData.UnmaskedAerosolBackscatterCoefficient*(770/828);
    HSRL.BSR828 = HSRL.Ba828./HSRL.Bm828+1;

    LidarData.OfflineCombinedTotalCounts = Counts.foff;
    LidarData.OfflineMolecularTotalCounts = Counts.foff_mol;
    [LidarData]=BackscatterRetrievalRayleighBrillouin070621(LidarData,WeatherData);
    HSRL.BSRf = LidarData.UnmaskedBackscatterRatio;

%     LidarData.OfflineCombinedTotalCounts = Counts.goff;
%     LidarData.OfflineMolecularTotalCounts = Counts.goff_mol;
%     [LidarData]=BackscatterRetrievalRayleighBrillouin070621(LidarData,WeatherData);
%     HSRL.BSRg = LidarData.UnmaskedBackscatterRatio;
elseif span_days(1)>=datetime(2021,3,12,'TimeZone','UTC')
    LidarData.Range = Range.rm;
    LidarData.Time = Time.ts;
    LidarData.OfflineCombinedAverageCounts = Counts.o2off;
    LidarData.OfflineMolecularAverageCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    [LidarData]=BackscatterRetrievalRayleighBrillouin0312(LidarData,WeatherData);
    HSRL.BSR = LidarData.BackscatterRatio;
    HSRL.Ba = LidarData.AerosolBackscatterCoefficient;
    HSRL.Bm = LidarData.MolecularBackscatterCoefficient;
    Counts.o2off_mol_corrected = LidarData.o2off_mol_corr;

    LidarData.OfflineCombinedAverageCounts = Counts.foff-Counts.foff_bg;
    LidarData.OfflineMolecularAverageCounts = Counts.foff_mol-Counts.foff_mol_bg;
    [LidarData]=BackscatterRetrievalRayleighBrillouin0312(LidarData,WeatherData);
    HSRL.BSRf = LidarData.BackscatterRatio;

    LidarData.OfflineCombinedAverageCounts = Counts.goff-Counts.foff_bg;
    LidarData.OfflineMolecularAverageCounts = Counts.goff_mol-Counts.foff_mol_bg;
    [LidarData]=BackscatterRetrievalRayleighBrillouin0312(LidarData,WeatherData);
    HSRL.BSRg = LidarData.BackscatterRatio;
elseif span_days(1)>datetime(2021,2,20,'TimeZone','UTC')
    LidarData.Range = Range.rm;
    LidarData.Time = Time.ts;
    LidarData.OfflineCombinedAverageCounts = Counts.o2off;
    LidarData.OfflineMolecularAverageCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    [LidarData]=BackscatterRetrievalRayleighBrillouin(LidarData,WeatherData);
    HSRL.BSR = LidarData.BackscatterRatio; 
    Counts.o2off_mol_corrected = LidarData.o2off_mol_corr;
elseif span_days(1)>datetime(2021,1,28,'TimeZone','UTC')
    load('Overlap1104.mat','overlapcorrection')
    overlapcorrection = interp1(Range.rm_raw_o2,overlapcorrection,Range.rm);
    LidarData.OfflineCombinedAverageCounts = Counts.o2off./overlapcorrection;
    LidarData.OfflineMolecularAverageCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    [LidarData]=BackscatterRetrieval_2_10_21(LidarData,WeatherData);
    HSRL.BSR = LidarData.BackscatterRatio;
else
    load('Overlap1006.mat','overlapcorrection')
    overlapcorrection = interp1(rm_raw_o2,overlapcorrection,rm);
    [HSRL.Bm,HSRL.Ba,~]= BackscatterRatioV4(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
    load('Overlap1104.mat','overlapcorrection')
    overlapcorrection = interp1(rm_raw_o2,overlapcorrection,rm);
    %[Bm,Ba,BR]= BackscatterRetrievalRayleigh(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
    LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
    LidarData.OfflineMolecularAverageCounts = o2off_mol;
    WeatherData.Temperature = T;
    WeatherData.Pressure = P;
    [LidarData]=BackscatterRetrieval(LidarData,WeatherData);
    HSRL.BSR = LidarData.BackscatterRatio;
end

%%
%===== Calculate Model Counts =====
%   HSRL.BSR = ones(size(Counts.o2on));
%   HSRL.Ba = zeros(size(HSRL.Ba));
[Model] = modelCounts(Counts,Model,HSRL,Time,Range,Spectrum);