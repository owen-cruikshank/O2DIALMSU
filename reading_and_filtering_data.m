clear all

%%
%====================
%==== Constants =====
%====================
g0 = 9.80665;       %[m/s^2] Gravitational acceleration 
M_air = 0.0289644;  %[kg/mol] Molar mass of Earth's air 
R = 8.3144598;      %[J/(mol*K)] Universal gas constant 

c = 2.99792458E8;           %[m/s] Speed of light 
kb = 1.38065E-23;           %[J/K][m^2 kg s-2 K-1] Boltzman's constant 
h = 6.626E-34;              %[Js] Planck's constant 
mo2 = 5.314E-26;            %[kg] Mass O2 molecule 
mWV = 2.9915e-26;           %[kg] Mass H2O molecule
m_air = 4.792E-26;          %[kg] Mass of air
q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)
%%

disp('Reading in files')
 
%Read data from MSU DIAL
span_days = datetime(2020,5,9,'TimeZone','UTC');%yyyy,mm,dd

path = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\o2DIAL_data';

[DataStructure2, Options] = loadMSUNETcdf(span_days,path);
%%

disp('Creating time and range vectors')
% ====================
% === Time Vectors ===
% ====================
%  t_end = ts_raw_o2_on(end)-ts_raw_o2_on(1); %[s] Ending time  
%  t_step = 60;                               %[s] Time step
%  ts = 0:t_step:t_end;                       %[s] Time vector
%  thr = ts / 60 / 60;                        %[hr] Time vector
% ts = Options.TimeGrid*60*60;
% thr = Options.TimeGrid;
% % Shorten time grid to length of data if full day of data is not recorded
% ts = ts(1:size(DataStructure2.MCS.Channel2.Data,1));
% thr = thr(1:size(DataStructure2.MCS.Channel2.Data,1));


ts = DataStructure2.MCS.Channel2.NewTimeGrid*60*60;
thr = DataStructure2.MCS.Channel2.NewTimeGrid;
t_step = ts(2)-ts(1);                         %[s] Time step
i_time = length(ts);                          %[none] length of time vector


%t_pulse = 1;                    %[microseconds] Pulse duration

% =====================
% === Range Vectors ===
% =====================
% Create range vector from length of data bins
nsPerBin = double(250);                             %[ns] Nanoseconds per bin
NBins = double(560);                                %[none] Number of range bins
rangeBin = (c * nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
rangeMin = -150; %[m]
rangeMin = -200; %[m]
%rangeMin = 0; %[m]
rm_raw_o2 = rangeMin:rangeBin:NBins(1)*rangeBin+rangeMin-rangeBin;    %[m] Create range vector
rm_raw_o2 = rm_raw_o2(:);                           %[m] Convert range vector to column vector
r_max = 6000;                                       %[m] Max range 
rm = rm_raw_o2(rm_raw_o2<=r_max & rm_raw_o2>0);     %[m] Shorten range vector to max range
rkm = rm./1000;                                     %[km] Range vector
i_range = length(rm);                               %[none] Size of range vector

%%
%==============================
%== MSU Weather station data ==
%==============================

% Use if there is no data from weather station
%weather_Temperature_interp = ones(size(ts))*16;
%weather_absPressure_interp = ones(size(ts))*1000;
[weather_Temperature_interp, weather_absPressure_interp, weather_VW_interp] = ORSLweather(span_days,ts,path);

%%

% O2 Counts
o2on_intp = DataStructure2.MCS.Channel2.Data';
o2off_intp = DataStructure2.MCS.Channel10.Data';
o2on_intp_mol = DataStructure2.MCS.Channel0.Data';
o2off_intp_mol = DataStructure2.MCS.Channel8.Data';


% ======================
% === O2 Backscatter ===
% ======================
% --- O2 Background Subtraction ---
% Online
bg_o2on = mean(o2on_intp(end-20:end,:));% Take mean of last data points
o2on_bgsub = o2on_intp - bg_o2on;       % Background subtracted
o2on_bgsub(o2on_bgsub < 0) = 0;         % Minimum of zero

bg_o2on_mol = mean(o2on_intp_mol(end-20:end,:));% Take mean of last data points
o2on_bgsub_mol = o2on_intp_mol - bg_o2on_mol;       % Background subtracted
o2on_bgsub_mol(o2on_bgsub_mol < 0) = 0;         % Minimum of zero


% Offline
bg_o2off = mean(o2off_intp(end-20:end,:));% Take mean of last data points
o2off_bgsub = o2off_intp - bg_o2off;      % Background subtracted
o2off_bgsub(o2off_bgsub < 0) = 0;         % Minimum of zero

bg_o2off_mol = mean(o2off_intp_mol(end-20:end,:));% Take mean of last data points
o2off_bgsub_mol = o2off_intp_mol - bg_o2off_mol;      % Background subtracted
o2off_bgsub_mol(o2off_bgsub_mol < 0) = 0;         % Minimum of zero


% Interpolating to shorter range vector
o2on_noise = interp2(ts,rm_raw_o2,o2on_bgsub,ts,rm);
o2on_noise = fillmissing(o2on_noise,'nearest',1);
o2on_noise = fillmissing(o2on_noise,'nearest',2);

o2on_noise_mol = interp2(ts,rm_raw_o2,o2on_bgsub_mol,ts,rm);
o2on_noise_mol = fillmissing(o2on_noise_mol,'nearest',1);
o2on_noise_mol = fillmissing(o2on_noise_mol,'nearest',2);

o2off_noise = interp2(ts,rm_raw_o2,o2off_bgsub,ts,rm);
o2off_noise = fillmissing(o2off_noise,'nearest',1);
o2off_noise = fillmissing(o2off_noise,'nearest',2);

o2off_noise_mol = interp2(ts,rm_raw_o2,o2off_bgsub_mol,ts,rm);
o2off_noise_mol = fillmissing(o2off_noise_mol,'nearest',1);
o2off_noise_mol = fillmissing(o2off_noise_mol,'nearest',2);


%%
% --- O2 Filtering  ---

% Minutes to average data over
t_avg = 30;                     %[min]

% Range oversample
oversample = 8;                 %[bins] Oversample must be even

% Moving average in time and range
% Replicate edges to prevent problems caused by zero padding
% Even kernel causes data to shift by half step -- interpolate back to original time & range vectors
k = ones(oversample,t_avg)./(oversample*t_avg);     % Kernel

% Online
o2on_noise_pad = padarray(o2on_noise,[oversample/2,t_avg/2],'replicate');
o2on_filt = filter2(k,o2on_noise_pad,'valid');
o2on = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt(1:end-1,1:end-1),ts,rm);
o2on = fillmissing(o2on,'nearest',1); % Fill in NaNs in dimension 1
o2on = fillmissing(o2on,'nearest',2); % Fill in NaNs in dimension 2

o2on_noise_pad_mol = padarray(o2on_noise_mol,[oversample/2,t_avg/2],'replicate');
o2on_filt_mol = filter2(k,o2on_noise_pad_mol,'valid');
%o2on_filt_mol = o2on_noise_mol;
%o2on_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt_mol,ts,rm);
o2on_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt_mol(1:end-1,1:end-1),ts,rm);
o2on_mol = fillmissing(o2on_mol,'nearest',1); % Fill in NaNs in dimension 1
o2on_mol = fillmissing(o2on_mol,'nearest',2); % Fill in NaNs in dimension 2

% Offline
o2off_noise_pad = padarray(o2off_noise,[oversample/2,t_avg/2],'replicate');
o2off_filt = filter2(k,o2off_noise_pad,'valid');
o2off = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt(1:end-1,1:end-1),ts,rm);
o2off = fillmissing(o2off,'nearest',1); % Fill in NaNs in dimension 1
o2off = fillmissing(o2off,'nearest',2); % Fill in NaNs in dimension 2

o2off_noise_pad_mol = padarray(o2off_noise_mol,[oversample/2,t_avg/2],'replicate');
o2off_filt_mol = filter2(k,o2off_noise_pad_mol,'valid');
o2off_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt_mol(1:end-1,1:end-1),ts,rm);
o2off_mol = fillmissing(o2off_mol,'nearest',1); % Fill in NaNs in dimension 1
o2off_mol = fillmissing(o2off_mol,'nearest',2); % Fill in NaNs in dimension 2