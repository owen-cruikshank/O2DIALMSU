function [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadBoulderdata(span_days,Options,Constant)

cd ../
pathMATLAB = [pwd '\Data\NCAR Boulder Data\MatlabPreload\'];
pathPython = [pwd '\Data\NCAR Boulder Data\Python\'];
pathSonde = [pwd '\Data\NCAR Boulder Data\'];
cd ../../../
homepath = pwd;
%cd( '.\OneDrive - Montana State University - Bozeman\Research\O2 DIAL\analysis')
cd( '.\OneDrive - Montana State University\Research\O2 DIAL\analysis')

for jj = 1:length(span_days)
    [Data] = loadNCARBoulderData(span_days(jj),pathMATLAB);
    [DataPython] = loadNCARBoulderDataPython(span_days(jj),pathPython);
    DataSonde = loadNCARBoulderDataSonde(span_days(jj),pathSonde);

    T_sgp_raw = {[]};
    rm_sgp = {[]};
    P_sgp_raw = {[]};
    datetime_sgp_cell = {[]};

    % T_sgp_raw = DataSonde.temperature;
    % rm_sgp = DataSonde.rm;
    % P_sgp_raw = DataSonde.pressure;
    % datetime_sgp_cell = DataSonde.date;
    % t_single_sgp = DataSonde.date_ts;
    if jj==1
    ts_raw_o2_on = Data.Lidar.Interp.O2OfflineComb.TimeStamp' * 60*60;
    ts_raw_o2_off = ts_raw_o2_on;

    bins = 49;
    o2on_raw = Data.Lidar.Interp.O2OnlineComb.Data' /bins;
    o2off_raw = Data.Lidar.Interp.O2OfflineComb.Data' /bins;
    O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data = Data.Lidar.Interp.O2OfflineMol.Data' /bins;
    O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data = Data.Lidar.Interp.O2OnlineMol.Data' /bins;

    %rm_raw = rm_raw' - 46.8-113;
    range_Backscatter_Ratio = DataPython.range;
    time_Backscatter_Ratio = DataPython.time';
    Backscatter_Ratio = DataPython.Backscatter_Ratio;

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

    O2Offline_Wavelength = Data.TimeSeries.Laser.O2Offline.WavelengthActual';
    O2Online_Wavelength = Data.TimeSeries.Laser.O2Online.WavelengthActual';

    Counts.NBins = Data.Lidar.Raw.O2OfflineMol.NBins';
    else
        ts_raw_o2_on = [ts_raw_o2_on Data.Lidar.Interp.O2OfflineComb.TimeStamp'*60*60+ts_raw_o2_on(end)];
        ts_raw_o2_off = ts_raw_o2_on;

        bins = 49;
        o2on_raw = [o2on_raw Data.Lidar.Interp.O2OnlineComb.Data'/bins];
        o2off_raw = [o2off_raw Data.Lidar.Interp.O2OfflineComb.Data' /bins];
        O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data = [O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data Data.Lidar.Interp.O2OfflineMol.Data'/bins];
        O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data = [O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data Data.Lidar.Interp.O2OnlineMol.Data'/bins];

        %rm_raw = rm_raw' - 46.8-113;
        range_Backscatter_Ratio = DataPython.range;
        time_Backscatter_Ratio = [time_Backscatter_Ratio DataPython.time'+time_Backscatter_Ratio(end)];
        Backscatter_Ratio = [Backscatter_Ratio DataPython.Backscatter_Ratio];

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

        O2Offline_Wavelength = [O2Offline_Wavelength Data.TimeSeries.Laser.O2Offline.WavelengthActual'];
        O2Online_Wavelength = [O2Online_Wavelength Data.TimeSeries.Laser.O2Online.WavelengthActual'];

        Counts.NBins = [Counts.NBins Data.Lidar.Raw.O2OfflineMol.NBins'];
        
    end
end

%%
disp('Conditioning files')




%convert time from single to double
% ts_raw_o2_on = double(time_O2_Offline_Backscatter_Channel_Raw_Data);   %[s] Raw time vector from data files 
% ts_raw_o2_off = double(time_O2_Online_Backscatter_Channel_Raw_Data);   %[s] Raw time vector from data files 

% o2on_raw = double(O2_Online_Backscatter_Channel_Raw_Data);      
% o2off_raw = double(O2_Offline_Backscatter_Channel_Raw_Data);   


disp('Creating time and range vectors')
% ====================
% === Time Vectors ===
% ====================
t_start=ts_raw_o2_on(1);
 t_end = ts_raw_o2_on(end);%-ts_raw_o2_on(1); %[s] Ending time  
 Time.t_step = 60;                               %[s] Time step
 Time.t_step = 60*Options.intTime;                               %[s] Time step
% ts = t_start:t_step:t_end;                       %[s] Time vector
Time.ts = 0:Time.t_step:t_end;
 thr = Time.ts / 60 / 60;                        %[hr] Time vector
 Time.i_time = length(Time.ts);                       %[none] length of time vector

%%
%Fill in nans

[y,m,d]=ymd(span_days(1));
date_ts_N = datetime(y,m,d,0,0,Time.ts,'TimeZone','UTC');
[hours,minute,sec] = hms(date_ts_N);

ts_D_NLogical = ~((hours == 0) & (minute == 0));

ts_D_N = Time.ts(ts_D_NLogical); 

Time.date_ts = date_ts_N;

index = 1;
for j = 1:length(date_ts_N)
    if index<=length(Time.date_ts)
        if isequal(datenum(date_ts_N(j)),datenum(Time.date_ts(index)))
            mask_nodata(:,j) = ones(133,1);
            index = index+1;
        else
            mask_nodata(:,j) = zeros(133,1);
        end
    else
        mask_nodata(:,j) = zeros(133,1);
    end
end


% =====================
% === Range Vectors ===
% =====================
%create range vector from length of data bins
Range.nsPerBin = double(250);                             %[ns] Convert data to double
%Range.NBins = double(560);                                %[none] Convert data to double
Range.NBins = double(490);  
Range.rangeBin = (Constant.c * Range.nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
%rangeBin = rangeBin*2;
Range.rm_raw_o2 = Range.rangeBin:Range.rangeBin:Range.NBins(1)*Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector

% rm_raw = 250e-9*(0:length(o2on_raw(:,1))-1)*Constant.c/2;
% rm_raw = rm_raw' - 46.8;

Range.r_max = 5000;                                       %[m] Max range 
Range.rm = Range.rm_raw_o2(Range.rm_raw_o2<=Range.r_max);                   %[m] Shorten range vector

%Range.rm = Range.rm(1:2:end);%integrate to new
%Range.rangeBin = Range.rangeBin*2;

Range.rkm = Range.rm./1000;                                     %[km] Range vector
Range.i_range = length(Range.rm);                               %[none] Size of range vector

%%
% Count dead time correction
nsPerBin = 250e-9;
summedBins = 14000/2;
%summedBins = 560;
% Deadtime for apd SPCM-780-12-FC
deadTime = 22e-9;
o2on_countRate = o2on_raw /(nsPerBin * summedBins);
o2on_correctionFactor = 1./(1 - deadTime*o2on_countRate);
o2on_raw_corrected = o2on_correctionFactor .* o2on_raw;
o2off_countRate = o2off_raw /(nsPerBin * summedBins);
o2off_correctionFactor = 1./(1 - deadTime*o2off_countRate);
o2off_raw_corrected = o2off_correctionFactor .* o2off_raw;

Counts.o2on_raw_corrected = o2on_raw_corrected;
Counts.o2off_raw_corrected = o2off_raw_corrected;

% O2 Counts
% Interpolating counts to match new time vector
% % % % o2on_intp = interp2(ts_raw_o2_on,rm_raw,o2on_raw_corrected,Time.ts,Range.rm_raw_o2,'nearest',nan);
% % % % o2off_intp = interp2(ts_raw_o2_off,rm_raw,o2off_raw_corrected,Time.ts,Range.rm_raw_o2,'nearest',nan);
% % % % o2off_intp_mol = interp2(ts_raw_o2_off,rm_raw,double(O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data),Time.ts,Range.rm_raw_o2,'nearest',nan);
% % % % o2on_intp_mol = interp2(ts_raw_o2_off,rm_raw,double(O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data),Time.ts,Range.rm_raw_o2,'nearest',nan);

%integrate to new time grid
[o2on_intp,o2onNBins] = intSum(o2on_raw_corrected,ts_raw_o2_on,Counts.NBins,Time.ts);
[o2off_intp,o2offNBins] = intSum(o2off_raw_corrected,ts_raw_o2_on,Counts.NBins,Time.ts);
[o2off_intp_mol,o2offNBins] = intSum(double(O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data),ts_raw_o2_on,Counts.NBins,Time.ts);
[o2on_intp_mol,o2offNBins] = intSum(double(O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data),ts_raw_o2_on,Counts.NBins,Time.ts);

%interpolate to new range
% o2on_intp = interp2(Time.ts,Range.rm_raw_o2-37.4741/8,o2on_intp,Time.ts,Range.rm_raw_o2,'linear',nan);
% o2off_intp = interp2(Time.ts,Range.rm_raw_o2,o2off_intp,Time.ts,Range.rm_raw_o2,'linear',nan);
% o2on_intp_mol = interp2(Time.ts,Range.rm_raw_o2,o2on_intp_mol,Time.ts,Range.rm_raw_o2,'linear',nan);
% o2off_intp_mol = interp2(Time.ts,Range.rm_raw_o2,o2off_intp_mol,Time.ts,Range.rm_raw_o2,'linear',nan);

%integrate to new range
% inc = 1;
% for ii = 1:2:length(Range.rm_raw_o2)
%     o2on_intp2(inc,:) = (o2on_intp(ii,:)+o2on_intp(ii+1,:))/2;
%     o2off_intp2(inc,:) = (o2off_intp(ii,:)+o2off_intp(ii+1,:))/2;
%     o2on_intp2_mol(inc,:) = (o2on_intp_mol(ii,:)+o2on_intp_mol(ii+1,:))/2;
%     o2off_intp2_mol(inc,:) = (o2off_intp_mol(ii,:)+o2off_intp_mol(ii+1,:))/2;
%     inc = inc+1;
% end
% o2on_intp = o2on_intp2;
% o2off_intp = o2off_intp2;
% o2on_intp_mol = o2on_intp2_mol;
% o2off_intp_mol = o2off_intp2_mol;

%Range.rm_raw_o2 = Range.rm_raw_o2(1:2:end);


%%
% ======================
% === O2 Backscatter ===
% ======================
Counts.o2on_raw = o2on_raw;
Counts.o2off_raw = o2off_raw;
% --- O2 Background Subtraction ---

endPoints = 150;
% Online
Counts.bg_o2on = mean(o2on_intp(end-endPoints:end,:),1,'omitnan');% Take mean of last data points
Counts.o2on_bgsub = o2on_intp - Counts.bg_o2on;       % Background subtracted
Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
Counts.o2on_noise = Counts.o2on_bgsub(1:Range.i_range,:);   % Shorten to length of range vector

% Offline
Counts.bg_o2off = mean(o2off_intp(end-endPoints:end,:),1,'omitnan');% Take mean of last data points
Counts.o2off_bgsub = o2off_intp - Counts.bg_o2off;      % Background subtracted
Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
Counts.o2off_noise = Counts.o2off_bgsub(1:Range.i_range,:);   % Shorten to length of range vector

% Online
Counts.bg_o2on_mol = mean(o2on_intp_mol(end-endPoints:end,:),1,'omitnan');% Take mean of last data points
Counts.o2on_bgsub_mol = o2on_intp_mol - Counts.bg_o2on_mol;       % Background subtracted
Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
Counts.o2on_noise_mol = Counts.o2on_bgsub_mol(1:Range.i_range,:);   % Shorten to length of range vector

% Offline
Counts.bg_o2off_mol = mean(o2off_intp_mol(end-endPoints:end,:),1,'omitnan');% Take mean of last data points
Counts.o2off_bgsub_mol = o2off_intp_mol - Counts.bg_o2off_mol;      % Background subtracted
Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero
Counts.o2off_noise_mol = Counts.o2off_bgsub_mol(1:Range.i_range,:);   % Shorten to length of range vector

%%
% load('AfterPulse.mat','pulseON','pulseOFF','pulseON_mol','pulseOFF_mol')
% Counts.o2on_noise(5:end,:) = Counts.o2on_noise(5:end,:)-pulseON(5:133,:)*-1;
% Counts.o2off_noise(5:end,:) = Counts.o2off_noise(5:end,:)-pulseOFF(5:133,:)*-1;
% Counts.o2on_noise_mol(5:end,:) = Counts.o2on_noise_mol(5:end,:)-pulseON_mol(5:133,:)*-1;
% Counts.o2off_noise_mol(5:end,:) = Counts.o2off_noise_mol(5:end,:)-pulseOFF_mol(5:133,:)*-1;
% 
% Counts.o2on_noise(Counts.o2on_noise<0)=0;
% Counts.o2off_noise(Counts.o2off_noise<0)=0;
% Counts.o2on_noise_mol(Counts.o2on_noise_mol<0)=0;
% Counts.o2off_noise_mol(Counts.o2off_noise_mol<0)=0;

%%
% --- O2 Filtering ---
% Moving average in time and range
% Replicate edges to prevent problems caused by zero padding
% Even kernel causes data to shift by half step -- interpolate back to original time & range vectors
k = ones(Options.oversample,Options.t_avg)./(Options.oversample*Options.t_avg);     % Kernel

% Online
% % o2on_noise_pad = padarray(o2on_noise,[oversample/2,t_avg/2],'replicate');
% % o2on_filt = filter2(k,o2on_noise_pad,'valid');
% % o2on = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt(1:end-1,1:end-1),ts,rm);
Counts.o2on = filter2(k,Counts.o2on_noise,'same');
Counts.o2on = fillmissing(Counts.o2on,'nearest',1); % Fill in NaNs in dimension 1
Counts.o2on = fillmissing(Counts.o2on,'nearest',2); % Fill in NaNs in dimension 2

% Offline
% % o2off_noise_pad = padarray(o2off_noise,[oversample/2,t_avg/2],'replicate');
% % o2off_filt = filter2(k,o2off_noise_pad,'valid');
% % o2off = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt(1:end-1,1:end-1),ts,rm);
Counts.o2off = filter2(k,Counts.o2off_noise,'same');
Counts.o2off = fillmissing(Counts.o2off,'nearest',1); % Fill in NaNs in dimension 1
Counts.o2off = fillmissing(Counts.o2off,'nearest',2); % Fill in NaNs in dimension 2

% % o2on_noise_pad_mol = padarray(o2on_noise_mol,[oversample/2,t_avg/2],'replicate');
% % o2on_filt_mol = filter2(k,o2on_noise_pad_mol,'valid');
% % o2on_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt_mol(1:end-1,1:end-1),ts,rm);
Counts.o2on_mol = filter2(k,Counts.o2on_noise_mol,'same');
Counts.o2on_mol = fillmissing(Counts.o2on_mol,'nearest',1); % Fill in NaNs in dimension 1
Counts.o2on_mol = fillmissing(Counts.o2on_mol,'nearest',2); % Fill in NaNs in dimension 2

% Offline
% % o2off_noise_pad_mol = padarray(o2off_noise_mol,[oversample/2,t_avg/2],'replicate');
% % o2off_filt_mol = filter2(k,o2off_noise_pad_mol,'valid');
% % o2off_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt_mol(1:end-1,1:end-1),ts,rm);
Counts.o2off_mol = filter2(k,Counts.o2off_noise_mol,'same');
Counts.o2off_mol = fillmissing(Counts.o2off_mol,'nearest',1); % Fill in NaNs in dimension 1
Counts.o2off_mol = fillmissing(Counts.o2off_mol,'nearest',2); % Fill in NaNs in dimension 2

%% 
%overlap correction between on and off
% load('OverlapONOFFBoulder.mat','OverlapONOFF')
% Counts.o2off = Counts.o2off.*OverlapONOFF;

%overlap correction
% load('overlapDiffCorr061621.mat','lnOonOoff')
% Counts.o2off = Counts.o2off.*smoothdata(lnOonOoff,1);

load('overlapCorr061821.mat','overlapCorr')
% % %Counts.o2off = Counts.o2off.*overlapCorr;
[~,constAlt]=min(abs(2000-Range.rm));

Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*overlapCorr(1:constAlt,:);
%%

%=====================
%= Backscatter ratio =
%=====================
%make backscatter ratio the same dimentions as data
%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio),double(Backscatter_Ratio),ts,rm); % interpolate backscatter ratio to range and time to match other data
%===interpolate to new time grid===
%%%%%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio),double(Backscatter_Ratio),Time.ts,Range.rm,'nearest',nan); % interpolate backscatter ratio to range and time to match other data
%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio)-3*(range_Backscatter_Ratio(2)-range_Backscatter_Ratio(1)),double(Backscatter_Ratio),ts,rm,'nearest',nan); % interpolate backscatter ratio to range and time to match other data

%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio)-1*(range_Backscatter_Ratio(2)-range_Backscatter_Ratio(1)),double(Backscatter_Ratio),ts,rm,'nearest',nan); % interpolate backscatter ratio to range and time to match other data

%===integrate to new time grid===
[BSRintSum,BSRbins] = intSum(double(Backscatter_Ratio),double(time_Backscatter_Ratio),Counts.NBins,Time.ts);
for ii = 1:Time.i_time
    HSRL.BSR(:,ii) = interp1(double(range_Backscatter_Ratio),BSRintSum(:,ii),Range.rm);
end
HSRL.BSR = fillmissing(HSRL.BSR,'nearest',1);
HSRL.BSR = fillmissing(HSRL.BSR,'nearest',2);

%BSR = BSR.*1.15;

%filter backscatter ratio
% % % BSR_noise_pad = padarray(BSR,[oversample/2,t_avg/2],'replicate');
% % % BSR_filt = filter2(k,BSR_noise_pad,'valid');
% % % BSR = interp2(ts-t_step/2,rm-rangeBin/2,BSR_filt(1:end-1,1:end-1),ts,rm);
% BSR = fillmissing(BSR,'nearest',1); % Fill in NaNs in dimension 1
% BSR = fillmissing(BSR,'nearest',2); % Fill in NaNs in dimension 2

% % %set end of BSR to 1
% % BSR = BSR./mean(BSR(end-10:end,:),1);

% BSR_unfiltered = BSR;
% BSR_filter = sgolayfilt(BSR, 1, 2*oversample+1);
% BSR = BSR_filter;

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
[abs_humidintSum,abs_humidbins] = intSum(double(Absolute_Humidity),double(time_Absolute_Humidity),Counts.NBins,Time.ts);
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


O2Online_Wavelength = interp1(time_Pressure,O2Online_Wavelength,Time.ts);
O2Offline_Wavelength = interp1(time_Pressure,O2Offline_Wavelength,Time.ts);

%lambda_online = 769.2330;                          %[nm] Cobleigh online wavelength
%lambda_offline = 769.3184;                         %[nm] Cobleigh offline wavelength
Spectrum.lambda_online = O2Online_Wavelength*10^9;          %[nm] Updated online wavelength
Spectrum.lambda_offline = O2Offline_Wavelength*10^9;        %[nm] Updated offline wavelength

Spectrum.lambda_online = 769.7958;          %[nm] Updated online wavelength
Spectrum.lambda_offline = 770.1085;        %[nm] Updated offline wavelength
Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber

%nu01 = nu_online(1000);                                   %[cm-1] Set center of scan to 
nu01 = Spectrum.nu_online(1);                                   %[cm-1] Set center of scan to 
nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
%nuBin = 0.00111;                                    %[cm-1] Scan increment
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
i_scan = length(nu_scan);                           %[none] length of scan vector

nuMin_off = Spectrum.nu_offline-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334;                                 %[cm-1] Scan upper bound
nu_scan_off = (nuMin_off:Spectrum.nuBin:nuMax_off);


lambda_scan = 10^7./nu_scan;                        %[nm](1 x lambda) Scan vector
%f_scan = nu_scan * c * 100;                         %[Hz]( x f) Scan vector

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.lambda_scan_3D_short = 10^7./Spectrum.nu_scan_3D_short;
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);         %[none] length of scan vector

Spectrum.del_nu = Spectrum.nu_scan_3D_short-Spectrum.nu_online;                %[1/cm] difference from center
%del_f = del_nu*100*c*10^-9;                         %[GHz] difference from center


[~,Spectrum.online_index] = min(abs(Spectrum.nu_online - Spectrum.nu_scan_3D_short),[],3);%finding index of online wavenumber
[~,Spectrum.offline_index] = min(abs(Spectrum.nu_offline - Spectrum.nu_scan_3D_short_off),[],3);%finding index of online wavenumber
tic
% % absorption = nan(i_range,i_time,1); %preallocate
% % for i = 1:i_range
% %     for j=1:i_time
% % absorption(i,j,:) = absorption_O2_770_model(T(i,j),P(i,j),nu_online(1,j),WV(i,j)); %[m-1] Function to calculate theoretical absorption
% %     end
% % end

Model.absorption(:,:,:) = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online(1,1),Model.WV); %[m-1] Function to calculate theoretical absorption

Model.transmission = exp(-cumtrapz(Range.rm,Model.absorption));
% % absorption_off = nan(i_range,i_time,1); %preallocate
% % for i = 1:i_range
% %     for j=1:i_time
% % absorption_off(i,j,:) = absorption_O2_770_model(T(i,j),P(i,j),nu_offline(1,j),WV(i,j)); %[m-1] Function to calculate theoretical absorption
% %     end
% % end

Model.absorption_off(:,:,:) = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline(1,1),Model.WV); %[m-1] Function to calculate theoretical absorption
toc
%absorption = absorption_O2_770_model_wavenumber(T,P,nu_online,WV); %[m-1] Funcrtion to calculate theoretical absorption
%[absorption_sonde,cross_section_sonde,lineshape_sonde] = absorption_O2_770_model(T_real,Patm_real,nu_online(1),0);
%[absorption_sonde,cross_section_sonde,lineshape_sonde] = absorption_O2_770_model_wavenumber(T_real,Patm_real,nu_online(1),0);

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

Data.MCS.Channel0.Data = nan(length(Range.rm_raw_o2),length(Time.ts));
Data.MCS.Channel2.Data = nan(length(Range.rm_raw_o2),length(Time.ts));
Data.MCS.Channel8.Data = nan(length(Range.rm_raw_o2),length(Time.ts));
Data.MCS.Channel10.Data = nan(length(Range.rm_raw_o2),length(Time.ts));

Data.Thermocouple.InsideCell.TimeStamp = [];
Data.Thermocouple.InsideCell.Temperature = [];
Data.Thermocouple.OutsideCell.Temperature = [];
Data.Thermocouple.TSOA.Temperature = [];
Data.Thermocouple.RoomTemp.Temperature = [];

%Data.Etalon.O2Etalon.TimeStamp = [];
%Data.Etalon.O2Etalon.TemperatureActual = [];
% Data.Laser.O2Online.TemperatureActual = [];
% Data.Laser.O2Online.TemperatureDesired = [];
% Data.Laser.O2Offline.TemperatureActual = [];
% Data.Laser.O2Offline.TemperatureDesired = [];
% Data.Laser.O2Online.TimeStamp = [];
% Data.Laser.O2Online.WaveDiff = [];
% Data.Laser.O2Offline.WaveDiff = [];
Data.Laser.TWSOA.TemperatureActual = [];
Data.Laser.TWSOA.TemperatureDesired = [];
Data.Laser.TWSOA.TimeStamp = [];

Data.UPS.all.TimeStamp = [];
Data.UPS.all.BatteryTimeLeft = [];
Data.UPS.all.HoursOnBattery = [];
Data.UPS.all.UPSTemperature = [];

%Reset Counts
Counts.NBins = o2onNBins;




