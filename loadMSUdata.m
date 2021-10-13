function [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data,Options] = loadMSUdata(span_days,Options,Constant)
%====================
%==== Reading Files =====
%====================
disp('Reading in files')

cd ../
Options.path = [pwd '\Data\MSU data\RSync\NetCDFOutput\']; %Path for instument netCDF data
Options.weatherPath = [pwd '\Data\']; %path for weather station data

cd ../../../
Options.sondepath = [pwd '\Box\Radiosondes\Data\All Data\']; %path for radiosonde data

cd( '.\OneDrive - Montana State University\Research\O2 DIAL\analysis') %back to main directory

Options.MPDname = 'MSU';
Options.BinTotal = 560;
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
Range.NBins = floor(560/Options.intRange); %number of range bins in vector
Range.rangeBin = (Constant.c * Range.nsPerBin(1)*10^-9)/2; %(m)range bin length

Range.rm_raw_o2 = -150:Range.rangeBin:Range.NBins(1)*Range.rangeBin-150-Range.rangeBin;    %[m] Create range vector
%Range.rm_raw_o2 = 0:Range.rangeBin:Range.NBins(1)*Range.rangeBin+0-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector
Range.r_max = 6000;                                       %[m] Max range 
Range.rm = Range.rm_raw_o2(Range.rm_raw_o2<=Range.r_max & Range.rm_raw_o2>0);     %[m] Shorten range vector to max range

%=== Integrate range vector ==
% Range.rm = Range.rm(1:2:end);%integrate to new
% Range.rangeBin = Range.rangeBin*2;

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
lapseRate = lapseRate / 1000;                       %[K/m] 

Model.T = Model.Ts + lapseRate .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
Model.P = Model.Ps .* (Model.Ts./Model.T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  

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
        if length(T_sonde_int{i})<Range.i_range % ==If sonde does not reach full lidar range
            disp('ran')
            Sonde.T_sonde(1:length(T_sonde_int{i}),i) = T_sonde_int{i};
            Sonde.T_sonde(length(T_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(T_sonde_int{i}),1);
            Sonde.P_sonde(1:length(T_sonde_int{i}),i) = P_sonde_int{i};
            Sonde.P_sonde(length(P_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(P_sonde_int{i}),1);
            Sonde.WV_sonde(1:length(T_sonde_int{i}),i) = WV_sonde_int{i};
            Sonde.WV_sonde(length(WV_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(WV_sonde_int{i}),1);
        else
            Sonde.T_sonde(:,i) = T_sonde_int{i}(1:Range.i_range);
            Sonde.P_sonde(:,i) = P_sonde_int{i}(1:Range.i_range);
            Sonde.WV_sonde(:,i) = WV_sonde_int{i}(1:Range.i_range);
        end
        %===interp sonde time
        [rm_sgp{i},IA,IC] = unique(rm_sgp{i});
        sonde_time(1:length(rm_sonde_int{i}),i) = interp1(rm_sgp{i},sondeStruc(i).time(IA),rm_sonde_int{i})';  
        if length(sonde_time) < Range.i_range
            sonde_time = [sonde_time; sonde_time(end).*ones(Range.i_range-length(sonde_time),1)];
        end
         %Find index of sonde in time vector
         for j = 1:Range.i_range
            [~, Sonde.sonde_ind(j,i)]=min(abs(sonde_datetime(i)+seconds(sonde_time(j,i))-Time.date_ts));
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
%%
%==== After pulse correction ====
% load('AfterPulse2.mat','PulseOn','PulseOff','PulseOnMol','PulseOffMol')
% Data.MCS.Channel2.Data = Data.MCS.Channel2.Data - PulseOn+1;
% Data.MCS.Channel10.Data = Data.MCS.Channel10.Data - PulseOff+1;
% Data.MCS.Channel0.Data = Data.MCS.Channel0.Data - PulseOnMol+1;
% Data.MCS.Channel8.Data = Data.MCS.Channel8.Data - PulseOffMol+1;


%=========================
%=== Calculate background
%=========================
Counts.bg_o2off = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
Counts.o2off_bgsub = Data.MCS.Channel10.Data - Counts.bg_o2off;       % Background subtracted
Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero

Counts.bg_o2on = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
Counts.o2on_bgsub = Data.MCS.Channel2.Data - Counts.bg_o2on;       % Background subtracted
Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero

Counts.bg_o2on_mol = mean(Data.MCS.Channel0.Data(end-20:end,:));% Take mean of last data points
Counts.o2on_bgsub_mol = Data.MCS.Channel0.Data - Counts.bg_o2on_mol;       % Background subtracted
Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero

Counts.bg_o2off_mol = mean(Data.MCS.Channel8.Data(end-20:end,:));% Take mean of last data points
Counts.o2off_bgsub_mol = Data.MCS.Channel8.Data - Counts.bg_o2off_mol;       % Background subtracted
Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero

%%

% ========integrate to new range
% % inc = 1;
% % for ii = 1:2:length(Range.rm_raw_o2)
% %     o2on_intp2(inc,:) = (Counts.o2on_bgsub(ii,:)+Counts.o2on_bgsub(ii+1,:))/2;
% %     o2off_intp2(inc,:) = (Counts.o2off_bgsub(ii,:)+Counts.o2off_bgsub(ii+1,:))/2;
% %     o2on_intp2_mol(inc,:) = (Counts.o2on_bgsub_mol(ii,:)+Counts.o2on_bgsub_mol(ii+1,:))/2;
% %     o2off_intp2_mol(inc,:) = (Counts.o2off_bgsub_mol(ii,:)+Counts.o2off_bgsub_mol(ii+1,:))/2;
% %     inc = inc+1;
% % end
% % Counts.o2on_bgsub = o2on_intp2;
% % Counts.o2off_bgsub = o2off_intp2;
% % Counts.o2on_bgsub_mol = o2on_intp2_mol;
% % Counts.o2off_bgsub_mol = o2off_intp2_mol;
% % 
% % Range.rm_raw_o2 = Range.rm_raw_o2(1:2:end);

%%

% ===Interpolating to shorter range vector===
Counts.o2on_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub,Time.ts,Range.rm);
Counts.o2on_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub_mol,Time.ts,Range.rm);
Counts.o2off_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub,Time.ts,Range.rm);
Counts.o2off_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub_mol,Time.ts,Range.rm);

% Counts.NBins = Data.MCS.Channel0.NBins*2;
Counts.NBins = Data.MCS.Channel0.NBins;

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
%lambda_online = 769.2330;                          %[nm] Cobleigh online wavelength
%lambda_offline = 769.3184;                         %[nm] Cobleigh offline wavelength

%lambda_online = interp1(Options.TimeGrid,Data.Laser.O2Online.WavelengthActual,Time.ts/60/60);
%lambda_offline = interp1(Options.TimeGrid,Data.Laser.O2Offline.WavelengthActual,Time.ts/60/60);

Spectrum.lambda_online = 769.7958 *ones(size(Time.ts));
Spectrum.lambda_offline = 770.1085 *ones(size(Time.ts));

Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber

nuMin = Spectrum.nu_online-0.334;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nuMin_off = Spectrum.nu_offline-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334;                                 %[cm-1] Scan upper bound
nu_scan_off = (nuMin_off:Spectrum.nuBin:nuMax_off);

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension

Spectrum.lambda_scan_3D_short = 10^7./Spectrum.nu_scan_3D_short;
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);         %[none] length of scan vector

Spectrum.del_nu = Spectrum.nu_scan_3D_short-Spectrum.nu_online;                %[1/cm] difference from center
Spectrum.del_lambda = Spectrum.lambda_scan_3D_short-Spectrum.lambda_online;

[~,Spectrum.online_index] = min(abs(Spectrum.nu_online - Spectrum.nu_scan_3D_short),[],3);%finding index of online wavenumber
[~,Spectrum.offline_index] = min(abs(Spectrum.nu_offline - Spectrum.nu_scan_3D_short_off),[],3);%finding index of online wavenumber

%===== Calculate Model absorption from Model T and P =======
Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption
Model.transmission = exp(-cumtrapz(Range.rm,Model.absorption));

%===== Calucation Model absorption for radiosondes =======
for i=1:numel(sonde_datetime) 
        if isdatetime(sonde_datetime(i)) %Check if there are any sondes
            Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(Sonde.sonde_ind(:,i)),Model.WV(:,Sonde.sonde_ind(:,i)))); %[m-1] Funcrtion to calculate theoretical absorption
            Sonde.trasmission_sonde{i} = exp(-cumtrapz(Range.rm,Sonde.absorption_sonde{i})); %O2 transmission
        else
            Sonde.absorption_sonde{i} = nan(Range.i_range,1);
            Sonde.trasmission_sonde{i} = nan(Range.i_range,1);
            Sonde.T_sonde = nan(Range.i_range,1);
            Sonde.P_sonde = nan(Range.i_range,1);
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
%====== Calucate any appy optimal filtering based on Poisson thinning ====
Counts = poissonThin(Counts);

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
elseif span_days(1)>=datetime(2021,7,6,'TimeZone','UTC')
    %addpath '\BSR retrieval'
    %load('Overlap0304.mat')
    %LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
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
%    Counts.o2off_mol_corrected = LidarData.o2off_mol_corr;
elseif span_days(1)>=datetime(2021,3,12,'TimeZone','UTC')
    %addpath '\BSR retrieval'
    %load('Overlap0304.mat')
    %LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
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
elseif span_days(1)>datetime(2021,2,20,'TimeZone','UTC')
    %addpath '\BSR retrieval'
    %load('Overlap0304.mat')
    %LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
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
    %[Bm,Ba,BR]= BackscatterRetrievalRayleigh(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
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
[Model] = modelCounts(Counts,Model,HSRL,Time,Range,Spectrum);