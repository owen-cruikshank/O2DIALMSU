function [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadBoulderdata2(span_days,Options,Constant)

disp('Reading in files')
%Read data from MSU DIAL

%%%%load('6_25_21.mat','Data','Options')
cd ../
Options.path = [pwd '\Data\NCAR Boulder Data\RawData\'];
pathPython = [pwd '\Data\NCAR Boulder Data\Python\'];
pathSonde = [pwd '\Data\NCAR Boulder Data\'];
cd ../../../
cd( '.\OneDrive - Montana State University\Research\O2 DIAL\analysis')

Options.MPDname = 'Boulder';
Options.BinTotal = 490;
[Data, Options] = loadMSUNETcdf(span_days,Options);
%%

for jj = 1:length(span_days)
    [DataPython] = loadNCARBoulderDataPython(span_days(jj),pathPython);
    DataSonde = loadNCARBoulderDataSonde(span_days(jj),pathSonde);

    T_sgp_raw = {[]};
    rm_sgp = {[]};
    P_sgp_raw = {[]};
    datetime_sgp_cell = {[]};

    if jj==1

    range_Backscatter_Ratio = DataPython.range;
    time_Backscatter_Ratio = DataPython.time';
    Backscatter_Ratio = DataPython.Backscatter_Ratio;
    Aerosol_Backscatter_Coefficient = DataPython.Aerosol_Backscatter_Coefficient;
    Backscatter_Ratio_mask = DataPython.Backscatter_Ratio_mask;
    HSRLMolecular_RayleighBrillioun = DataPython.HSRLMolecular_RayleighBrillioun;
    r_freqency = DataPython.r_freqency;

    Absolute_Humidity = DataPython.Absolute_Humidity;
    Absolute_Humidity_mask = DataPython.Absolute_Humidity_mask;
    time_Absolute_Humidity = DataPython.time';
    range_Absolute_Humidity = DataPython.range;

    time_Temperature = DataPython.time';
    range_Temperature = DataPython.range;
    Temperature = DataPython.Temperature_Model;
    Surface_Temperature_HSRL = DataPython.Surface_Temperature';

    time_Pressure = DataPython.time';
    range_Pressure = DataPython.range;
    Pressure = DataPython.Pressure_Model;
    Surface_Pressure_HSRL = DataPython.Surface_Pressure';

    else
        range_Backscatter_Ratio = DataPython.range;
        time_Backscatter_Ratio = [time_Backscatter_Ratio DataPython.time'+time_Backscatter_Ratio(end)];
        Backscatter_Ratio = [Backscatter_Ratio DataPython.Backscatter_Ratio];
        Aerosol_Backscatter_Coefficient = [Aerosol_Backscatter_Coefficient DataPython.Aerosol_Backscatter_Coefficient];
        Backscatter_Ratio_mask = [Backscatter_Ratio_mask DataPython.Backscatter_Ratio_mask];
        %HSRLMolecular_RayleighBrillioun = [HSRLMolecular_RayleighBrillioun DataPython.HSRLMolecular_RayleighBrillioun];
        HSRLMolecular_RayleighBrillioun = cat(3,HSRLMolecular_RayleighBrillioun,DataPython.HSRLMolecular_RayleighBrillioun);
        r_freqency = [r_freqency DataPython.r_freqency];

        Absolute_Humidity = [Absolute_Humidity DataPython.Absolute_Humidity];
        Absolute_Humidity_mask = [Absolute_Humidity_mask DataPython.Absolute_Humidity_mask];
        time_Absolute_Humidity = [time_Absolute_Humidity DataPython.time'+time_Absolute_Humidity(end)];
        range_Absolute_Humidity = DataPython.range;

        time_Temperature = [time_Temperature DataPython.time'+time_Temperature(end)];
        range_Temperature = DataPython.range;
        Temperature = [Temperature DataPython.Temperature_Model];
        Surface_Temperature_HSRL = [Surface_Temperature_HSRL DataPython.Surface_Temperature'];

        time_Pressure = [time_Pressure DataPython.time'+time_Pressure(end)];
        range_Pressure = DataPython.range;
        Pressure = [Pressure DataPython.Pressure_Model];
        Surface_Pressure_HSRL = [Surface_Pressure_HSRL DataPython.Surface_Pressure'];

        
    end
end

Counts.NBinsPython = ones(size(time_Temperature));

Spectrum.HSRLMolecular_RayleighBrillioun = HSRLMolecular_RayleighBrillioun;
Spectrum.r_freqency = r_freqency;
HSRL.Backscatter_Ratio_mask = Backscatter_Ratio_mask;


%%
disp('Creating time and range vectors')
% ====================
% === Time Vectors ===
% ====================
%Time.ts = Options.TimeGrid*60*60;
Time.ts = Data.MCS.Channel0.NewTimeGrid*60*60; 
Time.t_step = Time.ts(2)-Time.ts(1);
Time.i_time = length(Time.ts);
Time.date_ts = span_days(1) + seconds(Time.ts);
Time.date_ts.TimeZone = 'UTC';


% =====================
% === Range Vectors ===
% =====================
% Create range vector from length of data bins
Range.nsPerBin = 250; %[ns] bin length in nanosections
%Range.NBins = floor(560/Options.intRange); %number of range bins
Range.NBins = floor(490/Options.intRange); %number of range bins
Range.rangeBin = (Constant.c * Range.nsPerBin(1)*10^-9)/2; %range bin length


Range.rm_raw_o2 = 0:Range.rangeBin:Range.NBins(1)*Range.rangeBin+0-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector
Range.r_max = 6000;                                       %[m] Max range 
Range.rm = Range.rm_raw_o2(Range.rm_raw_o2<=Range.r_max & Range.rm_raw_o2>0);     %[m] Shorten range vector to max range


Range.rm = Range.rm(1:2:end);%integrate to new
Range.rangeBin = Range.rangeBin*2;

% % % Range.rm = Range.rm(1:4:end);%integrate to new
% % % Range.rangeBin = Range.rangeBin*4;

Range.i_range = length(Range.rm);                               %[none] Size of range vector
%%


disp('Calculating model')

% Calculating temperature and pressure model
            
Model.Ts = interp1(time_Temperature,Surface_Temperature_HSRL,Time.ts,'nearest',nan); %[K] Surface temperature             
Model.Ps = interp1(time_Pressure,Surface_Pressure_HSRL,Time.ts,'nearest',nan); %[atm] Absolute surface pressure

lapseRate = -6.5;                                   %[K/km] Guess adiabatic lapse rate  typically -6.5 up to 10km
lapseRate = lapseRate / 1000;                       %[K/m] 

%T = Ts + lapseRate .* rm;                           %[K] (1 x r) Temperature model as a function of r 

Model.T = interp2(double(time_Temperature),double(range_Temperature),double(Temperature),Time.ts,Range.rm,'nearest'); % Interpolate Temperature model to range and time to match other data
Model.T = fillmissing(Model.T,'nearest',1);
Model.T = fillmissing(Model.T,'nearest',2);
%%%%T(1,:) = Ts;


%P = Ps .* (Ts./T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  

Model.P = interp2(double(time_Pressure),double(range_Pressure),double(Pressure),Time.ts,Range.rm,'nearest'); % Interpolate pressure model to range and time to match other data
Model.P = fillmissing(Model.P,'nearest',1);
Model.P = fillmissing(Model.P,'nearest',2);
%%%P(1,:) = Ps;

%===============
%= Water Vapor =
%===============

% Profile
Absolute_Humidity = Absolute_Humidity.*(1-double(Absolute_Humidity_mask));
Absolute_Humidity(Absolute_Humidity<0) = 0;

%make water vapor the same dimentions as data
%===interpolate to new time grid===
%%%%abs_humid = interp2(double(time_Absolute_Humidity),double(range_Absolute_Humidity),double(Absolute_Humidity),Time.ts,Range.rm,'nearest',nan); % interpolate water vapor to range and time to match other data
%===integrate to new time grid===
[abs_humidintSum,abs_humidbins] = intSum(double(Absolute_Humidity),double(time_Absolute_Humidity),Counts.NBinsPython,Time.ts);
%shorten to new rm
%%%abs_humid = abs_humid(1:Range.i_range,:);
for ii = 1:Time.i_time
    abs_humid(:,ii) = interp1(double(range_Absolute_Humidity),abs_humidintSum(:,ii),Range.rm);
end
abs_humid = fillmissing(abs_humid,'nearest',1);
abs_humid = fillmissing(abs_humid,'nearest',2);



WV_intp = abs_humid / 1000 / Constant.mWV;        %[molecules/m^3] water vapor number density

% Fill in missing lower data with repeated values
WV_cut = 9;
WV_intp_fill = [NaN((WV_cut),Time.i_time); WV_intp(WV_cut+1:end,:)];
WV_intp_fill = fillmissing(WV_intp_fill,'nearest');

Model.WV = WV_intp_fill;                      %[molecules/m^3] water vapor number density

Model.WV = WV_intp;
%WV = zeros(size(WV));

%%
%-----BSR
%===integrate to new time grid===
[BSRintSum,BSRbins] = intSum(double(Backscatter_Ratio),double(time_Backscatter_Ratio),Counts.NBinsPython,Time.ts);
[Aerosol_Backscatter_CoefficientintSum,BSRbins] = intSum(double(Aerosol_Backscatter_Coefficient),double(time_Backscatter_Ratio),Counts.NBinsPython,Time.ts);
for ii = 1:Time.i_time
    HSRL.BSR(:,ii) = interp1(double(range_Backscatter_Ratio),BSRintSum(:,ii),Range.rm);
    HSRL.Ba(:,ii) = interp1(double(range_Backscatter_Ratio),Aerosol_Backscatter_CoefficientintSum(:,ii),Range.rm);
end
HSRL.BSR = fillmissing(HSRL.BSR,'nearest',1);
HSRL.BSR = fillmissing(HSRL.BSR,'nearest',2);

HSRL.Ba = fillmissing(HSRL.Ba,'nearest',1);
HSRL.Ba = fillmissing(HSRL.Ba,'nearest',2);

HSRL.Bm = HSRL.Ba./(HSRL.BSR-1);%'m^(-1) sr^(-1)'



%%
% load('AfterPulse2.mat','PulseOn','PulseOff','PulseOnMol','PulseOffMol')
% Data.MCS.Channel2.Data = Data.MCS.Channel2.Data - PulseOn+1;
% Data.MCS.Channel10.Data = Data.MCS.Channel10.Data - PulseOff+1;
% Data.MCS.Channel0.Data = Data.MCS.Channel0.Data - PulseOnMol+1;
% Data.MCS.Channel8.Data = Data.MCS.Channel8.Data - PulseOffMol+1;

% --- O2 Background Subtraction ---
% Online
% % Counts.bg_o2on = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2on_bgsub = Data.MCS.Channel2.Data - Counts.bg_o2on;       % Background subtracted
% % Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
% % 
% % Counts.bg_o2off = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2off_bgsub = Data.MCS.Channel10.Data - Counts.bg_o2off;       % Background subtracted
% % Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
% % 
% % Counts.bg_o2on_mol = mean(Data.MCS.Channel0.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2on_bgsub_mol = Data.MCS.Channel0.Data - Counts.bg_o2on_mol;       % Background subtracted
% % Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
% % 
% % Counts.bg_o2off_mol = mean(Data.MCS.Channel8.Data(end-20:end,:));% Take mean of last data points
% % Counts.o2off_bgsub_mol = Data.MCS.Channel8.Data - Counts.bg_o2off_mol;       % Background subtracted
% % Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero

% Counts.bg_o2on = mean(Data.MCS.Channel8.Data(end-20:end,:));% Take mean of last data points
% Counts.o2on_bgsub = Data.MCS.Channel8.Data - Counts.bg_o2on;       % Background subtracted
% Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2off = mean(Data.MCS.Channel0.Data(end-20:end,:));% Take mean of last data points
% Counts.o2off_bgsub = Data.MCS.Channel0.Data - Counts.bg_o2off;       % Background subtracted
% Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2on_mol = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
% Counts.o2on_bgsub_mol = Data.MCS.Channel10.Data - Counts.bg_o2on_mol;       % Background subtracted
% Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2off_mol = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
% Counts.o2off_bgsub_mol = Data.MCS.Channel2.Data - Counts.bg_o2off_mol;       % Background subtracted
% Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero

Counts.bg_o2on = mean(Data.MCS.Channel9.Data(end-20:end,:));% Take mean of last data points
Counts.o2on_bgsub = Data.MCS.Channel9.Data - Counts.bg_o2on;       % Background subtracted
Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero

Counts.bg_o2off = mean(Data.MCS.Channel1.Data(end-20:end,:));% Take mean of last data points
Counts.o2off_bgsub = Data.MCS.Channel1.Data - Counts.bg_o2off;       % Background subtracted
Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero

Counts.bg_o2on_mol = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
Counts.o2on_bgsub_mol = Data.MCS.Channel10.Data - Counts.bg_o2on_mol;       % Background subtracted
Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero

Counts.bg_o2off_mol = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
Counts.o2off_bgsub_mol = Data.MCS.Channel2.Data - Counts.bg_o2off_mol;       % Background subtracted
Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero


% 
% Counts.bg_o2on = mean(Data.MCS.Channel0.Data(end-20:end,:));% Take mean of last data points
% Counts.o2on_bgsub = Data.MCS.Channel0.Data - Counts.bg_o2on;       % Background subtracted
% Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2off = mean(Data.MCS.Channel8.Data(end-20:end,:));% Take mean of last data points
% Counts.o2off_bgsub = Data.MCS.Channel8.Data - Counts.bg_o2off;       % Background subtracted
% Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2on_mol = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
% Counts.o2on_bgsub_mol = Data.MCS.Channel2.Data - Counts.bg_o2on_mol;       % Background subtracted
% Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
% 
% Counts.bg_o2off_mol = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
% Counts.o2off_bgsub_mol = Data.MCS.Channel10.Data - Counts.bg_o2off_mol;       % Background subtracted
% Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero
% 

%%
%integrate to new range
inc = 1;
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

Range.rm_raw_o2 = Range.rm_raw_o2(1:2:end);


%%

% Interpolating to shorter range vector
Counts.o2on_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub,Time.ts,Range.rm);
Counts.o2on_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub_mol,Time.ts,Range.rm);
Counts.o2off_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub,Time.ts,Range.rm);
Counts.o2off_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub_mol,Time.ts,Range.rm);

Counts.NBins = Data.MCS.Channel0.NBins;



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
%-----Create Spectrum vectors----
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

lambda_scan = 10^7./nu_scan;                        %[nm](1 x lambda) Scan vector

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension

Spectrum.lambda_scan_3D_short = 10^7./Spectrum.nu_scan_3D_short;
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);         %[none] length of scan vector

Spectrum.del_nu = Spectrum.nu_scan_3D_short-Spectrum.nu_online;                %[1/cm] difference from center
Spectrum.del_lambda = Spectrum.lambda_scan_3D_short-Spectrum.lambda_online;

[~,Spectrum.online_index] = min(abs(Spectrum.nu_online - Spectrum.nu_scan_3D_short),[],3);%finding index of online wavenumber
[~,Spectrum.offline_index] = min(abs(Spectrum.nu_offline - Spectrum.nu_scan_3D_short_off),[],3);%finding index of online wavenumber

Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption
Model.transmission = exp(-cumtrapz(Range.rm,Model.absorption));
%absorption_off = absorption_O2_770_model(T,P,nu_offline);

% for i=1:numel(sonde_datetime)
%         if isdatetime(sonde_datetime(i))
%             Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(Sonde.sonde_ind(:,i)),Model.WV(:,Sonde.sonde_ind(:,i)))); %[m-1] Funcrtion to calculate theoretical absorption
%             Sonde.trasmission_sonde{i} = exp(-cumtrapz(Range.rm,Sonde.absorption_sonde{i})); %O2 transmission
%         else
%             Sonde.absorption_sonde{i} = nan(Range.i_range,1);
%             Sonde.trasmission_sonde{i} = nan(Range.i_range,1);
%             Sonde.T_sonde = nan(Range.i_range,1);
%             Sonde.P_sonde = nan(Range.i_range,1);
%         end
% end

%%


%%
%Count filtering
k = ones(Options.oversample,Options.t_avg)./(Options.oversample*Options.t_avg);     % Kernel

Counts.o2on = filter2(k,Counts.o2on_noise,'same');
% Counts.o2on = fillmissing(Counts.o2on,'nearest',1); % Fill in NaNs in dimension 1
% Counts.o2on = fillmissing(Counts.o2on,'nearest',2); % Fill in NaNs in dimension 2

Counts.o2on_mol = filter2(k,Counts.o2on_noise_mol,'same');
% Counts.o2on_mol = fillmissing(Counts.o2on_mol,'nearest',1); % Fill in NaNs in dimension 1
% Counts.o2on_mol = fillmissing(Counts.o2on_mol,'nearest',2); % Fill in NaNs in dimension 2

Counts.o2off = filter2(k,Counts.o2off_noise,'same');
% Counts.o2off = fillmissing(Counts.o2off,'nearest',1); % Fill in NaNs in dimension 1
% Counts.o2off = fillmissing(Counts.o2off,'nearest',2); % Fill in NaNs in dimension 2

Counts.o2off_mol = filter2(k,Counts.o2off_noise_mol,'same');
% Counts.o2off_mol = fillmissing(Counts.o2off_mol,'nearest',1); % Fill in NaNs in dimension 1
% Counts.o2off_mol = fillmissing(Counts.o2off_mol,'nearest',2); % Fill in NaNs in dimension 2

%%
% % load('overlapCorrMSU6_18_21.mat','overlapCorr')
% load('overlapCorrMSU6_21_21.mat','overlapCorr')
% %Counts.o2off = Counts.o2off.*overlapCorr;

% load('overlapCorrMSU6_24_21.mat','overlapCorr')
% load('overlapCorrMSU6_24_21_1000.mat','overlapCorr')
%  [~,constAlt]=min(abs(2000-Range.rm));
%  Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*overlapCorr(1:constAlt,:);

% load('BoulderlnOonOoff6_16to6_29_21.mat','lnOonOoff')
% lnOonOoff = smoothdata(lnOonOoff,2,'movmean',10);
%   [~,constAlt]=min(abs(500-Range.rm));
%   Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*lnOonOoff(1:constAlt,:);
%%
%Calculating absorption due to radiosonde measurements
sonde_time = [];
if ~isempty(sonde_time)
    for i=1:numel(sonde_time(1,:))
            if isdatetime(sonde_time(1,i))
                %absorption_sonde{i} = diag(absorption_O2_770_model(T_real(:,i),Patm_real(:,i),nu_online(data_col_real(:,i)),WV(:,data_col_real(:,i)))); %[m-1] Function to calculate theoretical absorption
                Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(1),Model.WV(:,Sonde.sonde_ind(:,i)))); %[m-1] Function to calculate theoretical absorption
            else
                Sonde.absorption_sonde{i} = nan(i_range,1);
                Sonde.T_sonde = nan(i_range,1);
                Sonde.P_sonde = nan(i_range,1);
            end
    end
else
      Sonde.absorption_sonde{1} = nan(Range.i_range,1);
      Sonde.T_sonde = nan(Range.i_range,1);
      Sonde.P_sonde = nan(Range.i_range,1);
      Sonde.sonde_ind = [];
      
%       absorption_sonde = {};
%       T_sonde = [];
%       P_sonde = [];
end
%%
Data.Thermocouple.TSOA.Temperature = nan(size(Data.Thermocouple.InsideCell.Temperature));
Data.Thermocouple.OutsideCell.Temperature = nan(size(Data.Thermocouple.InsideCell.Temperature));


%%
[Model] = modelCounts(Counts,Model,HSRL,Time,Range,Spectrum);