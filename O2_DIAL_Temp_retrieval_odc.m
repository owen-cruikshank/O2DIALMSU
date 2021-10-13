%Reading data, calculating pertabative absorption, and iterative
%temperature for O2 DIAL data

%Owen Cruikshank
%February 16, 2019

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
m_air = 4.792E-26;           %[kg] Mass of air
q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)
%%

disp('Reading in files')
%====================
%=== Reading Data ===
%====================
%addpath('D:\Owen\OneDrive - Montana State University - Bozeman\research s19\matlab model\O2 dial model\MPD program\')
%addpath('C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\matlab model\O2 dial model\MPD program')
%addpath('C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\ARM data analysis\ARM DOE SGP Data')
%addpath('D:\Owen\OneDrive - Montana State University - Bozeman\research s19\ARM data analysis\ARM DOE SGP Data')
%addpath('D:\Owen\OneDrive - Montana State University - Bozeman\research s19\ARM data analysis\ARM DOE SGP Data\sondes')
%datapath = 'ARM DOE SGP Data\ForOwen-20200623T183315Z-001\ForOwen\Soundings';
addpath('C:\Users\oencr\OneDrive - Montana State University\Research\O2 DIAL\Data\ARM DOE SGP Data')
addpath('C:\Users\d98b386\OneDrive - Montana State University\Research\O2 DIAL\Data\ARM DOE SGP Data')

date_begin = datetime(2020,8,29);
date_end   = datetime(2020,8,29);

date_begin = datetime(2019,4,20);
date_end   = datetime(2019,4,20);
span_days = date_begin:date_end;        % Days in set [datetime]

% % % 
analysispath = pwd;

for i = 1:length(span_days)
 date_i = span_days(i);
    d = num2str(yyyymmdd(date_i));
    filenameMPD = insertAfter('wv_dial05..Python.nc',10,d(3:end));
    %filenameMPD = insertAfter('mpd05..Python.nc',6,d(3:end));
    
    cd ..
    cd('Data\ARM DOE SGP Data')
    
    filenameSGP = dir(['sgpsondewnpnC1.b1.',d,'.*']);
    
    cd(analysispath)
    %filenameSGP = dir(['\ARM DOE SGP Data\sgpsondewnpnC1.b1.',d,'.*']);
    %%filenameSGPTWR = dir(['sgp1twrmrC1.c1.',d,'.*']);
        % === MPD Raw Data ===
    % Import Aerosol backscatter, Temperature, Pressure, Water Vapor
    [ba_raw{i},bam{i},r_ba{i},t_ba{i},T_surf_raw{i},T_raw{i},r_T{i},t_T{i},P_surf_raw{i},P_raw{i},r_P{i},...
        t_P{i},WV_surf{i},WVm_surf{i},WV_raw{i},WVm{i},r_WV{i},t_WV{i}] = getMPDdata(filenameMPD); 
    % Import online and offline O2 counts
    [o2on_raw{i},r_o2on{i},t_o2on{i},o2off_raw{i},r_o2off{i},...
        t_o2off{i}] = getMPDO2data(filenameMPD);
    
    % Add multiples of 24 hrs (86400 seconds) to each proceeding day so time
    % increases monotonically (instead of starting over each day).
    t_T{i} = t_T{i} + (i-1)*86400;
    t_WV{i} = t_WV{i} + (i-1)*86400;
    t_ba{i} = t_ba{i} + (i-1)*86400;
    t_o2on{i} = t_o2on{i} + (i-1)*86400;%#ok<*SAGROW>
    
    
    % === SGP Raw Radiosonde Data ===
    %jj = +(i-1)*length(filenameSGP);    % Accounts for i increasing in proceeding loops
    for j = 1:length(filenameSGP)
        [T_sgp_raw{j,i},P_sgp_raw{j,i},rm_sgp{j,i},t_single_sgp{j,i}] = getSGPdata(filenameSGP(j).name);
        dot_loc = strfind(filenameSGP(j).name,'.');
        time_num_sgp{j,i} = filenameSGP(j).name(dot_loc(3)+1:dot_loc(end)-1);
        time_char_sgp{j,i} = datenum(time_num_sgp{j,i},'HHMMSS');
        time_sgp{j,i} = datestr(time_char_sgp{j,i},'HH:MM:SS');
        datetime_sgp_cell{j,i} = date_i + time_sgp{j,i};
    end
    
    if isempty(filenameSGP)
        T_sgp_raw = {nan(1,318)};
        P_sgp_raw = {nan(1,318)};
        rm_sgp = {nan(1,319)};
        t_single_sgp = {nan(1,319)};
        datetime_sgp_cell = {nan(1)};
    end

    % === SGP Tower Data ====
   %%% [T_02m{i}, T_25m{i}, T_60m{i}, P_02m{i}, P_25m{i}, P_60m{i}, t_sgptwr{i}] = getSGPTWRdata(filenameSGPTWR.name);
    
    
    
    
    
    %select data file
file = ['ARM DOE SGP Data\wv_dial05.' d(3:end) '.Python.nc'];
%file = ['ARM DOE SGP Data\mpd05.' d(3:end) '.Python.nc'];
%file = ['mpd05.' d(3:end) '.Python.nc'];

%read in relevant data
time_O2_Online_Backscatter_Channel_Raw_Data{i}                             = ncread(file,'time_O2_Online_Backscatter_Channel_Raw_Data');                           %[s](1440x1){single} seconds since midnight
time_O2_Online_Backscatter_Channel_Raw_Data{i} = time_O2_Online_Backscatter_Channel_Raw_Data{i} + (i-1)*86400;
range_O2_Online_Backscatter_Channel_Raw_Data{i}                            = ncread(file,'range_O2_Online_Backscatter_Channel_Raw_Data');                          %[m](560x1){single}
O2_Online_Backscatter_Channel_Raw_Data{i}                                  = ncread(file,'O2_Online_Backscatter_Channel_Raw_Data');                                %[photons](560x14400){single}Unpolarization Online O2-DIAL Backscatter Returns on the Combined Detector'
%O2_Online_Backscatter_Channel_Raw_Data_variance                         = ncread(file,'O2_Online_Backscatter_Channel_Raw_Data_variance');
%O2_Online_Backscatter_Channel_Raw_Data_ProfileCount                     = ncread(file,'O2_Online_Backscatter_Channel_Raw_Data_ProfileCount');                   %[count](1440x1){int}Number of raw profiles integrated into this profile
O2Offline_Wavelength{i}                                                    = ncread(file,'O2Offline_Wavelength');                                                  %[](1440x1){double} Time resolved Oxygen Offline wavelength
time_O2_Offline_Backscatter_Channel_Raw_Data{i}                            = ncread(file,'time_O2_Offline_Backscatter_Channel_Raw_Data');                          %[s](1440x1)
time_O2_Offline_Backscatter_Channel_Raw_Data{i} = time_O2_Offline_Backscatter_Channel_Raw_Data{i} + (i-1)*86400;
%range_O2_Offline_Backscatter_Channel_Raw_Data{i}                           = ncread(file,'range_O2_Offline_Backscatter_Channel_Raw_Data');                         %[m](560x1)
O2_Offline_Backscatter_Channel_Raw_Data{i}                                 = ncread(file,'O2_Offline_Backscatter_Channel_Raw_Data');                               %[photons](560x1440) Offline O2-DIAL Backscatter Returns on the Combined Detector
%O2_Offline_Backscatter_Channel_Raw_Data_variance                        = ncread(file,'O2_Offline_Backscatter_Channel_Raw_Data_variance');
%O2_Offline_Backscatter_Channel_Raw_Data_ProfileCount                    = ncread(file,'O2_Offline_Backscatter_Channel_Raw_Data_ProfileCount');                  %[count](1440x1){int}Number of raw profiles integrated into this profile
O2Online_Wavelength{i}                                                     = ncread(file,'O2Online_Wavelength');                                                   %[](1440x1){double} Time resolved Oxygen Offline wavelength
time_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{i}         = ncread(file,'time_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data');       %[s](1440x1){single} seconds since midnight
time_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{i} = time_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{i} + (i-1)*86400;
%range_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{i}        = ncread(file,'range_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data');      %[m](560x1){single}
O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{i}              = ncread(file,'O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data');            %[photons](560x1440) Offline O2-DIAL Backscatter Returns on the Molecular Detector
%O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data_variance     = ncread(file,'O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data_variance');   %[](560x1440){single}
%O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data_ProfileCount = ncread(file,'O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data_ProfileCount');%[count](1140x1){int}
time_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{i}          = ncread(file,'time_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data');        %[s](1440x1){single} seconds since midnight
time_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{i} = time_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{i} + (i-1)*86400;
%range_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{i}         = ncread(file,'range_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data');       %[m](560x1){single}
O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{i}               = ncread(file,'O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data');             %[photons](560x1440) Online O2-DIAL Backscatter Returns on the Molecular Detector
%O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data_variance      = ncread(file,'O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data_variance');    %[](560x1440){single}
%O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data_ProfileCount  = ncread(file,'O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data_ProfileCount');%[count](1140x1){int}

%time_Surface_Temperature_HSRL                                           = ncread(file,'time_Surface_Temperature_HSRL');
%range_Surface_Temperature_HSRL                                          = ncread(file,'range_Surface_Temperature_HSRL');
Surface_Temperature_HSRL{i}                                                = ncread(file,'Surface_Temperature_HSRL');
%Surface_Temperature_HSRL_variance                                       = ncread(file,'Surface_Temperature_HSRL_variance');
%Surface_Temperature_HSRL_ProfileCount                                   = ncread(file,'Surface_Temperature_HSRL_ProfileCount');
%time_Surface_Pressure_HSRL                                              = ncread(file,'time_Surface_Pressure_HSRL');
%range_Surface_Pressure_HSRL                                             = ncread(file,'range_Surface_Pressure_HSRL');
Surface_Pressure_HSRL{i}                                                   = ncread(file,'Surface_Pressure_HSRL');
%Surface_Pressure_HSRL_variance                                          = ncread(file,'Surface_Pressure_HSRL_variance');
%Surface_Pressure_HSRL_ProfileCount                                      = ncread(file,'Surface_Pressure_HSRL_ProfileCount');

time_Temperature_HSRL{i}                                                   = ncread(file,'time_Temperature_HSRL');
time_Temperature_HSRL{i} = time_Temperature_HSRL{i} + (i-1)*86400;
range_Temperature_HSRL{i}                                                  = ncread(file,'range_Temperature_HSRL');
Temperature_HSRL{i}                                                        = ncread(file,'Temperature_HSRL');
%Temperature_HSRL_variance{i}                                               = ncread(file,'Temperature_HSRL_variance');
%Temperature_HSRL_ProfileCount{i}                                           = ncread(file,'Temperature_HSRL_ProfileCount');

time_Temperature{i}                                                        = ncread(file,'time_Temperature');
time_Temperature{i} = time_Temperature{i} + (i-1)*86400;
range_Temperature{i}                                                       = ncread(file,'range_Temperature');
Temperature{i}                                                             = ncread(file,'Temperature');
%Temperature_variance                                                    = ncread(file,'Temperature_variance');
%Temperature_ProfileCount                                                = ncread(file,'Temperature_ProfileCount');

time_Pressure{i}                                                           = ncread(file,'time_Pressure');
time_Pressure{i} = time_Pressure{i} + (i-1)*86400;
range_Pressure{i}                                                          = ncread(file,'range_Pressure');
Pressure{i}                                                                = ncread(file,'Pressure');
%Pressure_variance                                                       = ncread(file,'Pressure_variance');
%Pressure_ProfileCount                                                   = ncread(file,'Pressure_ProfileCount');

time_Backscatter_Ratio{i}                                                  = ncread(file,'time_Backscatter_Ratio');
time_Backscatter_Ratio{i} = time_Backscatter_Ratio{i} + (i-1)*86400;
range_Backscatter_Ratio{i}                                                 = ncread(file,'range_Backscatter_Ratio');
Backscatter_Ratio{i}                                                       = ncread(file,'Backscatter_Ratio');%description         = 'Ratio of combined to molecular backscatter'
%Backscatter_Ratio_variance                                              = ncread(file,'Backscatter_Ratio_variance');
Backscatter_Ratio_mask{i}                                                  = ncread(file,'Backscatter_Ratio_mask');
%Backscatter_Ratio_ProfileCount                                          = ncread(file,'Backscatter_Ratio_ProfileCount');

time_Aerosol_Backscatter_Coefficient{i}                                                  = ncread(file,'time_Aerosol_Backscatter_Coefficient');
time_Aerosol_Backscatter_Coefficient{i} = time_Aerosol_Backscatter_Coefficient{i} + (i-1)*86400;
range_Aerosol_Backscatter_Coefficient{i}                                                 = ncread(file,'range_Aerosol_Backscatter_Coefficient');
Aerosol_Backscatter_Coefficient{i}                                                       = ncread(file,'Aerosol_Backscatter_Coefficient');%description         = 'Ratio of combined to molecular backscatter'
%Backscatter_Ratio_variance                                              = ncread(file,'Backscatter_Ratio_variance');
Aerosol_Backscatter_Coefficient_mask{i}                                                  = ncread(file,'Aerosol_Backscatter_Coefficient_mask');
%Backscatter_Ratio_ProfileCount                                          = ncread(file,'Backscatter_Ratio_ProfileCount');

time_Molecular_Backscatter_Coefficient{i}                                                  = ncread(file,'time_Molecular_Backscatter_Coefficient');
time_Molecular_Backscatter_Coefficient{i} = time_Molecular_Backscatter_Coefficient{i} + (i-1)*86400;
range_Molecular_Backscatter_Coefficient{i}                                                 = ncread(file,'range_Molecular_Backscatter_Coefficient');
Molecular_Backscatter_Coefficient{i}                                                       = ncread(file,'Molecular_Backscatter_Coefficient');%description         = 'Ratio of combined to molecular backscatter'
%Backscatter_Ratio_variance                                              = ncread(file,'Backscatter_Ratio_variance');
%Molecular_Backscatter_Coefficient_mask{i}                                                  = ncread(file,'Molecular_Backscatter_Coefficient_mask');
%Backscatter_Ratio_ProfileCount                                          = ncread(file,'Backscatter_Ratio_ProfileCount');

time_Absolute_Humidity{i}                                                  = ncread(file,'time_Absolute_Humidity');
time_Absolute_Humidity{i} = time_Absolute_Humidity{i} + (i-1)*86400;
range_Absolute_Humidity{i}                                                 = ncread(file,'range_Absolute_Humidity');
Absolute_Humidity{i}                                                       = ncread(file,'Absolute_Humidity');%[g/m^(3)]
%Absolute_Humidity_variance                                              = ncread(file,'Absolute_Humidity_variance');
Absolute_Humidity_mask{i}                                                  = ncread(file,'Absolute_Humidity_mask');
%Absolute_Humidity_ProfileCount                                          = ncread(file,'Absolute_Humidity_ProfileCount');
    

[y,m,d]=ymd(span_days(i));

    
date_ts{i} = datetime(y,m,d,0,0,time_O2_Online_Backscatter_Channel_Raw_Data{i}-(i-1)*86400);
    
end
%%
%concatenating files


O2_Online_Backscatter_Channel_Raw_Data = double(cat(2,O2_Online_Backscatter_Channel_Raw_Data{:}));
O2Offline_Wavelength = double(cat(1,O2Offline_Wavelength{:}))';
O2Online_Wavelength = double(cat(1,O2Online_Wavelength{:}))';
O2_Offline_Backscatter_Channel_Raw_Data = double(cat(2,O2_Offline_Backscatter_Channel_Raw_Data{:}));
O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data  = double(cat(2,O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{:}));
time_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data  = double(cat(1,time_O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{:}))';
O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data = double(cat(2,O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data{:}));
Surface_Temperature_HSRL = double(cat(2,Surface_Temperature_HSRL{:}));
Surface_Pressure_HSRL = double(cat(2,Surface_Pressure_HSRL{:}));
Temperature_HSRL = double(cat(2,Temperature_HSRL{:}));
Temperature = double(cat(2,Temperature{:}));
Pressure = double(cat(2,Pressure{:}));
Backscatter_Ratio = double(cat(2,Backscatter_Ratio{:}));
Backscatter_Ratio_mask = double(cat(2,Backscatter_Ratio_mask{:}));
Molecular_Backscatter_Coefficient = double(cat(2,Molecular_Backscatter_Coefficient{:}));
%Molecular_Backscatter_Coefficient_mask = double(cat(2,Molecular_Backscatter_Coefficient_mask{:}));
Aerosol_Backscatter_Coefficient = double(cat(2,Aerosol_Backscatter_Coefficient{:}));
Aerosol_Backscatter_Coefficient_mask = double(cat(2,Aerosol_Backscatter_Coefficient_mask{:}));
Absolute_Humidity = double(cat(2,Absolute_Humidity{:}));
Absolute_Humidity_mask = double(cat(2,Absolute_Humidity_mask{:}));

%time
time_O2_Online_Backscatter_Channel_Raw_Data = double(cat(1,time_O2_Online_Backscatter_Channel_Raw_Data{:}))';
time_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data  = double(cat(1,time_O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data{:}))';
time_Temperature = double(cat(1,time_Temperature{:}))';
time_Temperature_HSRL = double(cat(1,time_Temperature_HSRL{:}))';
time_Pressure = double(cat(1,time_Pressure{:}))';
time_Backscatter_Ratio = double(cat(1,time_Backscatter_Ratio{:}))';
time_Molecular_Backscatter_Coefficient = double(cat(1,time_Molecular_Backscatter_Coefficient{:}))';
time_Aerosol_Backscatter_Coefficient = double(cat(1,time_Aerosol_Backscatter_Coefficient{:}))';
time_Absolute_Humidity = double(cat(1,time_Absolute_Humidity{:}))';
time_O2_Offline_Backscatter_Channel_Raw_Data = double(cat(1,time_O2_Offline_Backscatter_Channel_Raw_Data{:}))';
date_ts = cat(1,date_ts{:})';

%range
rm_raw = double(range_O2_Online_Backscatter_Channel_Raw_Data{1});
range_Temperature_HSRL = double(range_Temperature_HSRL{1});
range_Temperature = double(range_Temperature{1});
range_Pressure = double(range_Pressure{1});
range_Backscatter_Ratio = double(range_Backscatter_Ratio{1});
range_Molecular_Backscatter_Coefficient = double(range_Molecular_Backscatter_Coefficient{1});
range_Aerosol_Backscatter_Coefficient = double(range_Aerosol_Backscatter_Coefficient{1});
range_Absolute_Humidity = double(range_Absolute_Humidity{1});



%%

% === SGP Radiosonde data ===
% Convert contents of each cell to double
%if isnan(T_sgp_raw)
    
%else
T_sgp_raw = cellfun(@(x)double(x), T_sgp_raw,'UniformOutput',false);
P_sgp_raw = cellfun(@(x)double(x), P_sgp_raw,'UniformOutput',false);
rm_sgp = cellfun(@(x)double(x), rm_sgp,'UniformOutput',false);
%end





%%
disp('Conditioning files')

%convert time from single to double
ts_raw_o2_on = double(time_O2_Offline_Backscatter_Channel_Raw_Data);   %[s] Raw time vector from data files 
ts_raw_o2_off = double(time_O2_Online_Backscatter_Channel_Raw_Data);   %[s] Raw time vector from data files 

o2on_raw = double(O2_Online_Backscatter_Channel_Raw_Data);      
o2off_raw = double(O2_Offline_Backscatter_Channel_Raw_Data);   
% 
%%
% 
% 
% cd ../
% pathMATLAB = [pwd '\Data\NCAR Boulder Data\MatlabPreload\'];
% pathPython = [pwd '\Data\NCAR Boulder Data\Python\'];
% pathSonde = [pwd '\Data\NCAR Boulder Data\'];
% cd ../../../
% homepath = pwd;
% %cd( '.\OneDrive - Montana State University - Bozeman\Research\O2 DIAL\analysis')
% cd( '.\OneDrive - Montana State University\Research\O2 DIAL\analysis')
% 
% [Data] = loadNCARBoulderData(span_days,pathMATLAB);
% [DataPython] = loadNCARBoulderDataPython(span_days,pathPython);
% DataSonde = loadNCARBoulderDataSonde(span_days,pathSonde);
% 
% T_sgp_raw = {[]};
% rm_sgp = {[]};
% P_sgp_raw = {[]};
% datetime_sgp_cell = {[]};
% 
% % T_sgp_raw = DataSonde.temperature;
% % rm_sgp = DataSonde.rm;
% % P_sgp_raw = DataSonde.pressure;
% % datetime_sgp_cell = DataSonde.date;
% % t_single_sgp = DataSonde.date_ts;
% 
% 
% 
% ts_raw_o2_on = Data.Lidar.Interp.O2OfflineComb.TimeStamp' * 60*60;
% ts_raw_o2_off = ts_raw_o2_on;
% 
% bins = 49;
% 
% o2on_raw = Data.Lidar.Interp.O2OnlineComb.Data' /bins;
% o2off_raw = Data.Lidar.Interp.O2OfflineComb.Data' /bins;
% O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data = Data.Lidar.Interp.O2OfflineMol.Data' /bins;
% O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data = Data.Lidar.Interp.O2OnlineMol.Data' /bins;
% rm_raw = 250e-9*(0:length(o2on_raw(:,1))-1)*c/2;
% rm_raw = rm_raw' - 46.8;
% 
% %rm_raw = rm_raw' - 46.8-113;
% 
% 
% 
% range_Backscatter_Ratio = DataPython.range;
% time_Backscatter_Ratio = DataPython.time;
% Backscatter_Ratio = DataPython.Backscatter_Ratio;
% 
% Absolute_Humidity = DataPython.Absolute_Humidity;
% Absolute_Humidity_mask = DataPython.Absolute_Humidity_mask;
% time_Absolute_Humidity = DataPython.time;
% range_Absolute_Humidity = DataPython.range;
% 
% time_Temperature = DataPython.time;
% range_Temperature = DataPython.range;
% Temperature = DataPython.Temperature_Model;
% Surface_Temperature_HSRL = DataPython.Surface_Temperature;
% 
% time_Pressure = DataPython.time;
% range_Pressure = DataPython.range;
% Pressure = DataPython.Pressure_Model;
% Surface_Pressure_HSRL = DataPython.Surface_Pressure;
% 
% O2Offline_Wavelength = Data.TimeSeries.Laser.O2Offline.WavelengthActual;
% O2Online_Wavelength = Data.TimeSeries.Laser.O2Online.WavelengthActual;


%%

disp('Creating time and range vectors')
% ====================
% === Time Vectors ===
% ====================
t_start=ts_raw_o2_on(1);
 t_end = ts_raw_o2_on(end);%-ts_raw_o2_on(1); %[s] Ending time  
 t_step = 60;                               %[s] Time step
% ts = t_start:t_step:t_end;                       %[s] Time vector
ts = 0:t_step:86340;
 thr = ts / 60 / 60;                        %[hr] Time vector
 i_time = length(ts);                       %[none] length of time vector

%%
%Fill in nans

[y,m,d]=ymd(span_days(1));
date_ts_N = datetime(y,m,d,0,0,ts);
[hours,minute,sec] = hms(date_ts_N);

ts_D_NLogical = ~((hours == 0) & (minute == 0));

ts_D_N = ts(ts_D_NLogical); 

date_ts = date_ts_N;

index = 1;
for j = 1:length(date_ts_N)
    if index<=length(date_ts)
        if isequal(datenum(date_ts_N(j)),datenum(date_ts(index)))
            mask_nodata(:,j) = ones(133,1);
            index = index+1;
        else
            mask_nodata(:,j) = zeros(133,1);
        end
    else
        mask_nodata(:,j) = zeros(133,1);
    end
end

%%


% Minutes to average data over
t_avg = 30;                     %[min]
%t_avg = 15;                     %[min]

% Range oversample
oversample = 8;                 %[bins] Oversample must be even
oversample = 4;                 %[bins] Oversample must be even

%oversample = 9;
%t_avg = 31;
%oversample = 1;
%t_avg = 5;
%t_pulse = 1;                    %[microseconds] Pulse duration

% =====================
% === Range Vectors ===
% =====================
%create range vector from length of data bins
nsPerBin = double(250);                             %[ns] Convert data to double
NBins = double(560);                                %[none] Convert data to double
rangeBin = (c * nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
%rangeBin = rangeBin*2;
rm_raw_o2 = rangeBin:rangeBin:NBins(1)*rangeBin;    %[m] Create range vector
rm_raw_o2 = rm_raw_o2(:);                           %[m] Convert range vector to column vector
r_max = 5000;                                       %[m] Max range 
rm = rm_raw_o2(rm_raw_o2<=r_max);                   %[m] Shorten range vector
rkm = rm./1000;                                     %[km] Range vector
i_range = length(rm);                               %[none] Size of range vector


%%
% SGP Radiosonde Data
% Can't interpolate with array with repeating values, so use Kevin's sonde preparation function
for i = 1:numel(T_sgp_raw)
    if ~isempty(T_sgp_raw{i}) & ~isnan(T_sgp_raw{i})
        % Subtract first range value (site elevation) from whole vector
        rm_sgp{i} = rm_sgp{i} - rm_sgp{i}(1);
        % Collect radiosonde surface measurements
        T_sgp_surf(i) = T_sgp_raw{i}(1);
        P_sgp_surf(i) = P_sgp_raw{i}(1);
        % Convert datetimes from cells to vector
        datetime_sgp(i) = datetime_sgp_cell{i};
        % Custom interpolation function
        [T_sgp_int{i},P_sgp_int{i},rm_sgp_int{i}] = interp_sonde(T_sgp_raw{i},P_sgp_raw{i},rm_sgp{i},rangeBin);
        T_sgp(:,i) = T_sgp_int{i}(1:length(rm));
        P_sgp(:,i) = P_sgp_int{i}(1:length(rm));
        
        %interp sonde time
        [rm_sgp{i},IA,IC] = unique(rm_sgp{i});
        
        sonde_time(1:length(rm),i) = interp1(rm_sgp{i},t_single_sgp{i}(IA),rm)';
        %sonde_time = sonde_time';
        
                 %Find index of sonde in time vector
         for j = 1:length((rm))
            [~, data_col_real(j,i)]=min(abs(sonde_time(j,i)-date_ts_N));
         end
    else
        T_sgp = [];
        P_sgp = [];
        datetime_sgp=datetime([],[],[]);
        sonde_time = [];
        data_col_real = [];
        
    end
end
T_real = T_sgp;
Patm_real = P_sgp;
datetime_real = datetime_sgp;
i_time_real = length(datetime_real);

% Testing to see if starting at 60 m will help.... ************************************************************************
%%%T_surf = T_60m;
%%%P_surf = P_60m;

% Find time lapsed since beginning datetime
time_lapsed = datetime_real - date_begin;
[hrs_real,mins_real,sec_real] = hms(time_lapsed);

% Indices where "real" data is available
%data_col_real = hrs_real*60 + mins_real + 1;

%%%%%[~,data_col_real] = min(abs(ts' - (hrs_real*60*60+mins_real*60)));


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

%%
% O2 Counts
% Interpolating counts to match new time vector
o2on_intp = interp2(ts_raw_o2_on,rm_raw,o2on_raw_corrected,ts,rm_raw_o2,'nearest',nan);
o2off_intp = interp2(ts_raw_o2_off,rm_raw,o2off_raw_corrected,ts,rm_raw_o2,'nearest',nan);
o2off_intp_mol = interp2(ts_raw_o2_off,rm_raw,double(O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data),ts,rm_raw_o2,'nearest',nan);
o2on_intp_mol = interp2(ts_raw_o2_off,rm_raw,double(O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data),ts,rm_raw_o2,'nearest',nan);



% ======================
% === O2 Backscatter ===
% ======================
% --- O2 Background Subtraction ---
% Online
bg_o2on = mean(o2on_intp(end-150:end,:),1,'omitnan');% Take mean of last data points
o2on_bgsub = o2on_intp - bg_o2on;       % Background subtracted
o2on_bgsub(o2on_bgsub < 0) = 0;         % Minimum of zero
o2on_noise = o2on_bgsub(1:i_range,:);   % Shorten to length of range vector

% Offline
bg_o2off = mean(o2off_intp(end-150:end,:),1,'omitnan');% Take mean of last data points
o2off_bgsub = o2off_intp - bg_o2off;      % Background subtracted
o2off_bgsub(o2off_bgsub < 0) = 0;         % Minimum of zero
o2off_noise = o2off_bgsub(1:i_range,:);   % Shorten to length of range vector

% Online
bg_o2on_mol = mean(o2on_intp_mol(end-150:end,:),1,'omitnan');% Take mean of last data points
o2on_bgsub_mol = o2on_intp_mol - bg_o2on_mol;       % Background subtracted
o2on_bgsub_mol(o2on_bgsub_mol < 0) = 0;         % Minimum of zero
o2on_noise_mol = o2on_bgsub_mol(1:i_range,:);   % Shorten to length of range vector

% Offline
bg_o2off_mol = mean(o2off_intp_mol(end-150:end,:),1,'omitnan');% Take mean of last data points
o2off_bgsub_mol = o2off_intp_mol - bg_o2off_mol;      % Background subtracted
o2off_bgsub_mol(o2off_bgsub_mol < 0) = 0;         % Minimum of zero
o2off_noise_mol = o2off_bgsub_mol(1:i_range,:);   % Shorten to length of range vector

 %%
% % %count deconvolutino
% % pulseLength = 150; %m
% % %N_pulse = o2on_bgsub;
% % %N_pulse(1:9,:)=0;
% % 
% % 
% % k = ones(oversample,t_avg)./(oversample*t_avg);     % Kernel
% % o2on_bgsub2 = filter2(k,o2on_bgsub,'same');
% % N_pulse = o2on_bgsub2;
% % o2on_bgsub2(1:5,:)=0;
% % N_pulse(1:5,:)=0;
% % [o2on_intp_decon,errorEst] = countDeconvolution(N_pulse,o2on_bgsub2,rm_raw_o2,pulseLength);
% % 
% % %o2on_intp2(1:9,1:length(o2on_intp_decon(1,:)))=-1*o2on_intp_decon(39*8+1:40*8+1,:);
% % offset =0;
% % o2on_bgsub3(1:5,1:length(o2on_intp_decon(1,:)))=-1*o2on_intp_decon(131*4+1+offset:132*4+offset+1,:);
% % [o2on_intp_decon,errorEst] = countDeconvolution(N_pulse,o2on_bgsub3,rm_raw_o2,pulseLength);
% % 
% % 
% % o2off_bgsub2 = filter2(k,o2off_bgsub,'same');
% % N_pulse = o2off_bgsub2;
% % o2off_bgsub2(1:5,:)=0;
% % N_pulse(1:5,:)=0;
% % [o2off_intp_decon,errorEst] = countDeconvolution(N_pulse,o2off_bgsub2,rm_raw_o2,pulseLength);
% % 
% % %o2off_intp2(1:9,1:length(o2off_intp_decon(1,:)))=-1*o2off_intp_decon(39*8+1:40*8+1,:);
% % offset =-4;
% % o2off_bgsub3(1:5,1:length(o2off_intp_decon(1,:)))=-1*o2off_intp_decon(131*4+1+offset:132*4+offset+1,:);
% % [o2off_intp_decon,errorEst] = countDeconvolution(N_pulse,o2off_bgsub3,rm_raw_o2,pulseLength);
% % 
% % o2on = real(o2on_intp_decon(1:i_range,:));   % Shorten to length of range vector
% % o2on(o2on<0)=0;
% % o2off = real(o2off_intp_decon(1:i_range,:));   % Shorten to length of range vector
% % o2off(o2off<0)=0;
% % 
% % for i =1:i_time
% % o2on(:,i) = sgolayfilt(o2on(:,i),0,5);
% % o2off(:,i) = sgolayfilt(o2off(:,i),0,5);
% % end
% % 
% % figure()
% % p_point = 330;
% % subplot(2,1,1)
% % plot(rm_raw_o2,o2on_bgsub2(:,p_point),rm_raw_o2,o2on_intp_decon(:,p_point))
% % hold on
% % plot(rm_raw_o2(131*4+1+offset:132*4+offset+1,:),o2on_intp_decon(131*4+1+offset:132*4+offset+1,p_point),'*')
% % legend('Counts','counts decon')
% % subplot(2,1,2)
% % plot(rm_raw_o2,o2off_bgsub2(:,p_point),rm_raw_o2,o2off_intp_decon(:,p_point))
% % hold on
% % plot(rm_raw_o2(131*4+1+offset:132*4+offset+1,:),o2off_intp_decon(131*4+1+offset:132*4+offset+1,p_point),'*')
% % legend('Counts','counts decon')

% figure()
% %imagesc(log(o2on_bgsub2(:,200:400).*rm_raw_o2.^2))
% imagesc(o2on_bgsub2(:,200:400).*rm_raw_o2.^2)
% figure()
% imagesc(log(abs(o2on_intp_decon(:,200:400).*rm_raw_o2.^2)))


%%
%O2 summing
%%o2on = o2on_noise


%%

% --- O2 Filtering ---
% Moving average in time and range
% Replicate edges to prevent problems caused by zero padding
% Even kernel causes data to shift by half step -- interpolate back to original time & range vectors
k = ones(oversample,t_avg)./(oversample*t_avg);     % Kernel

% Online
% % o2on_noise_pad = padarray(o2on_noise,[oversample/2,t_avg/2],'replicate');
% % o2on_filt = filter2(k,o2on_noise_pad,'valid');
% % o2on = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt(1:end-1,1:end-1),ts,rm);
o2on = filter2(k,o2on_noise,'same');
o2on = fillmissing(o2on,'nearest',1); % Fill in NaNs in dimension 1
o2on = fillmissing(o2on,'nearest',2); % Fill in NaNs in dimension 2

% Offline
% % o2off_noise_pad = padarray(o2off_noise,[oversample/2,t_avg/2],'replicate');
% % o2off_filt = filter2(k,o2off_noise_pad,'valid');
% % o2off = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt(1:end-1,1:end-1),ts,rm);
o2off = filter2(k,o2off_noise,'same');
o2off = fillmissing(o2off,'nearest',1); % Fill in NaNs in dimension 1
o2off = fillmissing(o2off,'nearest',2); % Fill in NaNs in dimension 2

% % o2on_noise_pad_mol = padarray(o2on_noise_mol,[oversample/2,t_avg/2],'replicate');
% % o2on_filt_mol = filter2(k,o2on_noise_pad_mol,'valid');
% % o2on_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt_mol(1:end-1,1:end-1),ts,rm);
o2on_mol = filter2(k,o2on_noise_mol,'same');
o2on_mol = fillmissing(o2on_mol,'nearest',1); % Fill in NaNs in dimension 1
o2on_mol = fillmissing(o2on_mol,'nearest',2); % Fill in NaNs in dimension 2

% Offline
% % o2off_noise_pad_mol = padarray(o2off_noise_mol,[oversample/2,t_avg/2],'replicate');
% % o2off_filt_mol = filter2(k,o2off_noise_pad_mol,'valid');
% % o2off_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt_mol(1:end-1,1:end-1),ts,rm);
o2off_mol = filter2(k,o2off_noise_mol,'same');
o2off_mol = fillmissing(o2off_mol,'nearest',1); % Fill in NaNs in dimension 1
o2off_mol = fillmissing(o2off_mol,'nearest',2); % Fill in NaNs in dimension 2

% % % o2on=o2on_noise;
% % % o2off=o2off_noise;
% % % o2on_mol=o2on_noise_mol;
% % % o2off_mol=o2off_noise_mol;

% order = 2;
% o2off_sg = sgolayfilt(o2off_noise, order, t_avg+1,[],2);
% o2off = sgolayfilt(o2off_sg, order, oversample+1,[],1);
% o2off(o2off<0)=0;
% 
% o2off_sg_mol = sgolayfilt(o2off_noise_mol, order, t_avg+1,[],2);
% o2off_mol = sgolayfilt(o2off_sg_mol, order, oversample+1,[],1);
% o2off_mol(o2off_mol<0)=0;
% 
% o2on_sg = sgolayfilt(o2on_noise, order, t_avg+1,[],2);
% o2on = sgolayfilt(o2on_sg, order, oversample+1,[],1);
% o2on(o2on<0)=0;
% 
% o2on_sg_mol = sgolayfilt(o2on_noise_mol, order, t_avg+1,[],2);
% o2on_mol = sgolayfilt(o2on_sg_mol, order, oversample+1,[],1);
% o2on_mol(o2on_mol<0)=0;

%%
% % % % % %2D sg filtering
% % % % % order = 0;
% % % % % %sgfilt=sgsf_2d(-4:4,-15:15,order,order,0);
% % % % % sgfilt=sgsf_2d(-1*floor(t_avg/2):floor(t_avg/2),-1*floor(oversample/2):floor(oversample/2),order,order,0);
% % % % % %sgfilt=sgsf_2d(-10:10,-3:3,order,order,0);
% % % % % 
% % % % % 
% % % % % o2off = filter2(sgfilt,o2off_noise,'same');
% % % % % o2off(o2off<0)=0;
% % % % % 
% % % % % o2on = filter2(sgfilt,o2on_noise,'same');
% % % % % o2on(o2on<0)=0;
% % % % % 
% % % % % o2off_mol = filter2(sgfilt,o2off_noise_mol,'same');
% % % % % o2off_mol(o2off_mol<0)=0;
% % % % % 
% % % % % o2on_mol = filter2(sgfilt,o2on_noise_mol,'same');
% % % % % o2on_mol(o2on_mol<0)=0;

%%
%overlap adition
% % load('overlapFile4_20_19.mat')
% % o2on_mol = o2on_mol .* overlapFile;
% % o2on = o2on + o2on_mol;

%%

%=====================
%= Backscatter ratio =
%=====================
%make backscatter ratio the same dimentions as data
%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio),double(Backscatter_Ratio),ts,rm); % interpolate backscatter ratio to range and time to match other data
BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio),double(Backscatter_Ratio),ts,rm,'nearest',nan); % interpolate backscatter ratio to range and time to match other data
%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio)-3*(range_Backscatter_Ratio(2)-range_Backscatter_Ratio(1)),double(Backscatter_Ratio),ts,rm,'nearest',nan); % interpolate backscatter ratio to range and time to match other data

%BSR = interp2(double(time_Backscatter_Ratio),double(range_Backscatter_Ratio)-1*(range_Backscatter_Ratio(2)-range_Backscatter_Ratio(1)),double(Backscatter_Ratio),ts,rm,'nearest',nan); % interpolate backscatter ratio to range and time to match other data
BSR = fillmissing(BSR,'nearest',1);
BSR = fillmissing(BSR,'nearest',2);

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
abs_humid = interp2(double(time_Absolute_Humidity),double(range_Absolute_Humidity),double(Absolute_Humidity),ts,rm,'nearest',nan); % interpolate water vapor to range and time to match other data
abs_humid = fillmissing(abs_humid,'nearest',1);
abs_humid = fillmissing(abs_humid,'nearest',2);



WV_intp = abs_humid / 1000 / mWV;        %[molecules/m^3] water vapor number density

% Fill in missing lower data with repeated values
WV_cut = 9;
WV_intp_fill = [NaN((WV_cut),i_time); WV_intp(WV_cut+1:end,:)];
WV_intp_fill = fillmissing(WV_intp_fill,'nearest');

WV = WV_intp_fill;                      %[molecules/m^3] water vapor number density

WV = WV_intp*5;
%WV = zeros(size(WV));

%%

%======================
%= Cloud and SNR mask =
%======================
disp('Calculating masks')

cloud_p_point = 18;
%cloud_p_point = 3750/60;
SNR_threshold = 4;
SNR_threshold = 3.5;
SD_threshold = 6*10^9;
SD_threshold = 10^8;
SD_threshold = 15*10^8;
SD_threshold = 10*10^7; % ARM data

SNR_threshold = 4;
SD_threshold = 5*10^9; %Boulder data

SD_threshold = 3*10^8;
[SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2(o2on,o2off,rm,ts,cloud_p_point,SNR_threshold,SD_threshold,oversample,t_avg);




%%
disp('Calculating model')

% Calculating temperature and pressure model
            
Ts = interp1(time_Temperature,Surface_Temperature_HSRL,ts,'nearest',nan); %[K] Surface temperature             
Ps = interp1(time_Pressure,Surface_Pressure_HSRL,ts,'nearest',nan); %[atm] Absolute surface pressure

lapseRate = -6.5;                                   %[K/km] Guess adiabatic lapse rate  typically -6.5 up to 10km
lapseRate = lapseRate / 1000;                       %[K/m] 

%T = Ts + lapseRate .* rm;                           %[K] (1 x r) Temperature model as a function of r 

T = interp2(double(time_Temperature),double(range_Temperature),double(Temperature),ts,rm,'nearest'); % Interpolate Temperature model to range and time to match other data
T = fillmissing(T,'nearest',1);
T = fillmissing(T,'nearest',2);
%%%%T(1,:) = Ts;


%P = Ps .* (Ts./T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  

P = interp2(double(time_Pressure),double(range_Pressure),double(Pressure),ts,rm,'nearest'); % Interpolate pressure model to range and time to match other data
P = fillmissing(P,'nearest',1);
P = fillmissing(P,'nearest',2);
%%%P(1,:) = Ps;


O2Online_Wavelength = interp1(time_Pressure,O2Online_Wavelength,ts);
O2Offline_Wavelength = interp1(time_Pressure,O2Offline_Wavelength,ts);

%lambda_online = 769.2330;                          %[nm] Cobleigh online wavelength
%lambda_offline = 769.3184;                         %[nm] Cobleigh offline wavelength
lambda_online = O2Online_Wavelength*10^9;          %[nm] Updated online wavelength
lambda_offline = O2Offline_Wavelength*10^9;        %[nm] Updated offline wavelength

lambda_online = 769.7958;          %[nm] Updated online wavelength
lambda_offline = 770.1085;        %[nm] Updated offline wavelength
nu_online = 10^7./lambda_online;                    %[cm-1] Online wavenumber
nu_offline = 10^7./lambda_offline;                  %[cm-1] Offline wavenumber

%nu01 = nu_online(1000);                                   %[cm-1] Set center of scan to 
nu01 = nu_online(1);                                   %[cm-1] Set center of scan to 
nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
nuBin = 0.00222;                                    %[cm-1] Scan increment
%nuBin = 0.00111;                                    %[cm-1] Scan increment
nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
i_scan = length(nu_scan);                           %[none] length of scan vector

lambda_scan = 10^7./nu_scan;                        %[nm](1 x lambda) Scan vector
f_scan = nu_scan * c * 100;                         %[Hz]( x f) Scan vector

nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
lambda_scan_3D_short = 10^7./nu_scan_3D_short;
i_scan_3D_short = length(nu_scan_3D_short);         %[none] length of scan vector

del_nu = nu_scan_3D_short-nu_online;                %[1/cm] difference from center
del_f = del_nu*100*c*10^-9;                         %[GHz] difference from center


[~,online_index] = min(abs(nu_online - nu_scan_3D_short),[],3);%finding index of online wavenumber
[~,offline_index] = min(abs(nu_offline - nu_scan_3D_short),[],3);%finding index of online wavenumber

tic
% % absorption = nan(i_range,i_time,1); %preallocate
% % for i = 1:i_range
% %     for j=1:i_time
% % absorption(i,j,:) = absorption_O2_770_model(T(i,j),P(i,j),nu_online(1,j),WV(i,j)); %[m-1] Function to calculate theoretical absorption
% %     end
% % end

absorption(:,:,:) = absorption_O2_770_model(T,P,nu_online(1,1),WV); %[m-1] Function to calculate theoretical absorption

% % absorption_off = nan(i_range,i_time,1); %preallocate
% % for i = 1:i_range
% %     for j=1:i_time
% % absorption_off(i,j,:) = absorption_O2_770_model(T(i,j),P(i,j),nu_offline(1,j),WV(i,j)); %[m-1] Function to calculate theoretical absorption
% %     end
% % end

absorption_off(:,:,:) = absorption_O2_770_model(T,P,nu_offline(1,1),WV); %[m-1] Function to calculate theoretical absorption
toc
%absorption = absorption_O2_770_model_wavenumber(T,P,nu_online,WV); %[m-1] Funcrtion to calculate theoretical absorption
%[absorption_sonde,cross_section_sonde,lineshape_sonde] = absorption_O2_770_model(T_real,Patm_real,nu_online(1),0);
%[absorption_sonde,cross_section_sonde,lineshape_sonde] = absorption_O2_770_model_wavenumber(T_real,Patm_real,nu_online(1),0);

%Calculating absorption due to radiosonde measurements
if ~isempty(sonde_time)
    for i=1:numel(sonde_time(1,:))
            if isdatetime(sonde_time(1,i))
                %absorption_sonde{i} = diag(absorption_O2_770_model(T_real(:,i),Patm_real(:,i),nu_online(data_col_real(:,i)),WV(:,data_col_real(:,i)))); %[m-1] Function to calculate theoretical absorption
                absorption_sonde{i} = diag(absorption_O2_770_model(T_real(:,i),Patm_real(:,i),nu_online(1),WV(:,data_col_real(:,i)))); %[m-1] Function to calculate theoretical absorption
            else
                absorption_sonde{i} = nan(i_range,1);
                T_real = nan(i_range,1);
                Patm_real = nan(i_range,1);
            end
    end
else
      absorption_sonde{1} = nan(i_range,1);
      T_real = nan(i_range,1);
      Patm_real = nan(i_range,1);
      
%       absorption_sonde = {};
%       T_sonde = [];
%       P_sonde = [];
end



%%

%==========================
%= Pertabative absorption =
%==========================

disp('Calculating absorption')


% === Zeroth Order ===
ind_r_lo = 1:i_range-oversample;                                            % High range vector
ind_r_hi = 1+oversample:i_range;                                            % Low range vector

%ind_r_lo = 1:i_range-1;                                            % High range vector
%ind_r_hi = 2:i_range;                                            % Low range vector
ln_o2 = log((o2on(ind_r_lo,:).*o2off(ind_r_hi,:))./(o2on(ind_r_hi,:).*o2off(ind_r_lo,:))); % Natural log of counts

% Zeroth order term
alpha_0_raw = ln_o2./2./(rangeBin*oversample);                              %[1/m] 
%alpha_0_raw = ln_o2./2./(rangeBin);                              %[1/m] 



% --- Filter alpha_0 ---
% % Savitzky-Golay filtering
% alpha_0_sgfilt = sgolayfilt(alpha_0_raw, 1, 2*oversample+1);

% Moving average
% % alpha_0_pad = padarray(alpha_0_raw,[oversample/2,t_avg/2],'replicate');
% % alpha_0_filt = filter2(k,alpha_0_pad,'valid');
% % alpha_0 = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% alpha_0 = fillmissing(alpha_0,'nearest',1);                                 % Fill in NaNs in dimension 1
% alpha_0 = fillmissing(alpha_0,'nearest',2);                                 % Fill in NaNs in dimension 2



alpha_0 = interp2(ts,rm(ind_r_lo),alpha_0_raw,ts,rm);


% %     disp('SG derivative')
% %     int_der = -log(o2on)+log(o2off);
% %     tic
% % %    [~,g] = sgolay(2,oversample);
% %     [~,g] = sgolay(1,oversample);
% %     parfor j=1:i_time
% %         alpha_0(:,j) = conv(int_der(:,j), factorial(1)/(-rangeBin)^1 * g(:,2), 'same')/2;
% %     end
% %     toc



    
    alpha_0(alpha_0==Inf | alpha_0==-Inf)=0;
    alpha_0 = fillmissing(alpha_0,'nearest');


% % ------ fitting alpha zero ---------
% % alpha_0_fit = zeros(i_range,i_time);
% %  for i = 1:i_time
% %      alpha_0_fit_obj = fit(rm(26:106,1),alpha_0(26:106,i),'poly2');
% %      alpha_0_fit(:,i) = alpha_0_fit_obj(rm);
% %  end

%first order
fsr_O2 = 157.9;%[GHz] etalon free spectral range
finesse_O2 = 15.43;%etalon finesse
%filter_bw = 13;%[nm]bandwidth of filter
%filter_bw2 = 1;%[nm]bandwidth of filter
%NBF_blocking = 10^-6;%NBF out of band blocking
%NBF_blocking2 = 10^-4;%NBS out of band blocking
%NBF_transmission = 0.9;
%NBF_transmission2 = 0.7;

% ======================================
% === Normalized Etalon Transmission ===
% ======================================
% --- O2 channel ---
fsr_O2_nm = lambda_online.^2 * fsr_O2 / c;                          %[nm] Free spectral range converted to nm
delta_phase = 2 * pi * lambda_online.^2 ./ fsr_O2_nm ./ lambda_scan_3D_short; %[radians] Phase difference between orders
delta_phase = delta_phase - 2 * pi * lambda_online ./ fsr_O2_nm;     %[radians] Shifting phase so that zero is at the online wavelength
F_O2 = 1 ./ (sin(pi/2/finesse_O2)^2);                               %[none] Coefficient of finesse
T_etalon = 1 ./ (1 + F_O2 * sin(delta_phase/2).^2);                 %[none] Etalon transmission as a function of wavelength
%%
% Calculating purtabative absorption
altitude = .314;
[alpha_1, alpha_2] = pertAbsorption(alpha_0, T_etalon, T, P, rm,ts,rkm,m_air, nu_online, nu_scan_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,t_avg,c,kb,altitude);

% % ----- smoothing alpha 2 -----
% alpha_2s = zeros(i_range,i_time);
% for i = 1:i_time
%     alpha_2s(:,i) = 1.15.*smooth(alpha_2(:,i), 4*oversample+1);
% end

% === Total alpha ===
alpha_total_raw = alpha_0 + alpha_1 + alpha_2;
%alpha_total_raw = alpha_0 ;

% Force total alpha to its modeled surface value
[~,cut] = min(abs(rm-500));             % Index where rm is closest to chosen value
[~,cut] = min(abs(rm-0));             % Index where rm is closest to chosen value
%[~,cut] = min(abs(rm-1000));             % Index where rm is closest to chosen value
alpha_total = [absorption(1,:); NaN((cut - 1),i_time); alpha_total_raw(cut:end,:)];
alpha_total = fillmissing(alpha_total,'linear');
alpha_total = alpha_total(2:end,:);     % Remove surface value since rm starts at del_r

% Moving average
% % % alpha_total_pad = padarray(alpha_total,[oversample/2,t_avg/2],'replicate');
% % % alpha_total_filt = filter2(k,alpha_total_pad,'valid');
% % % alpha_total = interp2(ts-t_step/2,rm-rangeBin/2,alpha_total_filt(1:end-1,1:end-1),ts,rm);
% % % alpha_total = fillmissing(alpha_total,'nearest',1); % Fill in NaNs in dimension 1
% % % alpha_total = fillmissing(alpha_total,'nearest',2); % Fill in NaNs in dimension 2


% % ----- smoothing alpha total -----
% % % for i = 1:i_time
% % %     alpha_total(:,i) = smooth(alpha_total(:,i), 4*oversample+1);
% % % end

%%
%apply SNR mask
alpha_0m = alpha_0 .* SNRm;
alpha_0m(alpha_0m == 0) = NaN;                  % Replace mask with NaNs

%alpha_1 = alpha_1 .* SNRm;
%alpha_1(alpha_1 == 0) = NaN;                    % Replace mask with NaNs


%alpha_total_raw = alpha_total_raw .* SNRm;
%alpha_total_raw(alpha_total_raw == 0) = NaN;    % Replace mask with NaNs


alpha_totalm = alpha_total .* SNRm .* cloud_SDm_above;
%alpha_totalm(alpha_totalm <= 0) = NaN;          % Replace mask with NaNs

%%alpha_totalm = alpha_total;
% for i = 1:i_time
%     alpha_totalm(:,i) = smooth(alpha_totalm(:,i), 4*oversample+1);
% end

% for i = 1:i_time
%     alpha_totalm(:,i) = smooth(alpha_totalm(:,i), 4*oversample+1);
% end

% alpha_totalm = alpha_totalm .* SNRm .* cloud_SDm_above;
% alpha_totalm(alpha_totalm <= 0) = NaN;          % Replace mask with NaNs

alpha_0m = alpha_0m .* SNRm .* cloud_SDm_above;
alpha_0m(alpha_0m <= 0) = NaN;          % Replace mask with NaNs

%%
disp('Temp retrieval function')
%WV = 0; % set water vapor concentration to zero
tic
[T_final_test,L_fit_sm_test,Ts_fit,Patm_final,mean_lapse_rate,exclusion] =  temperatureRetrieval(T,ts,rm,P,WV,nu_online,alpha_totalm,SNRm,cloud_SDm_above);
toc

%%

k = ones(8,30)./(8*30);     % Kernel

k = ones(14,10)./(14*10);     % Kernel

Smoothing = repmat(gaussmf(linspace(-1,1,8)',[0.5,0]),1,30);
Smoothing = repmat(gaussmf(linspace(-1,1,14)',[0.5,0]),1,10);
Smoothing = Smoothing./sum(sum(Smoothing));

k=Smoothing;


%T_final_test_pad = padarray(T_final_test,[14/2,10/2],'replicate');
T_final_test_pad = padarray(T_final_test,[14/2,10/2],'replicate');
T_final_test_filt = filter2(k,T_final_test_pad,'valid');
T_final_tests = interp2(ts-t_step/2,rm-rangeBin/2,T_final_test_filt(1:end-1,1:end-1),ts,rm);
T_final_tests = fillmissing(T_final_tests,'nearest',1); % Fill in NaNs in dimension 1
T_final_tests = fillmissing(T_final_tests,'nearest',2); % Fill in NaNs in dimension 2




T_finalm = T_final_tests .* SNRm .* cloud_SDm_above .* mask_nodata;
T_finalm(T_finalm <= 0) = NaN;

[~,cut] = min(abs(rm-500));             % Index where rm is closest to chosen value
T_finalm = [ NaN((cut - 1),i_time); T_finalm(cut:end,:)];


 %%
% % % % %savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\6_26_20\';
% % % % savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\10_14_20\';
% % % % savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\10_13_20\';
% % % % savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\10_21_20\';
% % % % savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\10_24_20\';
% % % % savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\12_14_20\';
% savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\12_14_20\newCB\';
% savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\12_15_20\';
%savePath = 'C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\1_25_21\new\';
savePath = 'C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\1_27_21\';

 savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\2_18_21\';
 savePath = 'C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\3_29_21\newOwenProgram\';
 savePath = 'C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\3_29_21\newOwenProgramNoDelta\';
 savePath = 'C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\3_29_21\newnewOwenProgramNoDelta\';
 savePath = 'C:\Users\d98b386\OneDrive - Montana State University\research s19\Reports\3_29_21\newnewnewOwenProgramNoDelta\';
 
 savedata=0;
 if savedata ==1
for i = 1:length(span_days(1))
    [y,m,d]=ymd(span_days(i));
    y=num2str(y);
    m=num2str(m);
    d=num2str(d);
    if length(m)<2 
        m=['0' m];
    end
    if length(d)<2 
        d=['0' d];
    end
    date = [y m d];
    
%     save([savePath date],'alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
%         ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
%         ,'absorption_f','mask_nodata'...
%         ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
%         ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
%         ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')
    
        save([savePath date],'alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
        ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
        ,'mask_nodata'...
        ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
        ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
        ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')
end
 end
%
% load('D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\6_23_20\2019\20190417.mat','alpha_totalm'...
%         ,'alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','rkm','o2on','o2off','o2on_mol','o2off_mol'...
%         ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
%         ,'absorption_f','mask_nodata'...
%         ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
%         ,'o2on_noise','o2on_noise_mol','o2off_noise','o2off_noise_mol'...
%         ,'date_ts_N','data_col_real','date_ts','WV','Ts','exclusion','T_final_test','T_real','absorption_sonde')

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
tickAngle = 0;

% plot_time = 338;                         %[min];
% [~,p_point] = min(abs(plot_time-ts/60)); % Find closest value to 338min for comparison to other program
% p_point = 338;3.7

t_index =2;    % Choose which measurement to compare. Max value here is number of "real" measurements
  p_point = data_col_real(:,t_index);

plot_time = 8;                         %[hr];
%[~,p_point] = min(abs(plot_time-ts/60/60)); % Find closest value to 338min for comparison to other program
%p_point = p_point*ones(i_range,1);



figure(728590)
plot(datenum(date_ts_N),permute(L_fit_sm_test(1,:,:),[3 2 1]))
hold on
plot(ts([p_point(1) p_point(end)])/60/60,[-10 -4])
hold off
grid on

xlim([datenum(date_ts_N(1)) datenum(date_ts_N(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
xtickangle(tickAngle)
ylim([-11e-3 -3e-3]) 
title('Fitted lapse rate')
xlabel('Time (UTC)')
ylabel('Lapse rate (K/km)')

figure(728591)
%plot(ts/60/60,Ts - permute(Ts_fit(1,:,:),[3 2 1]))
plot(ts/60/60,T(1,:) - permute(Ts_fit(1,:,:),[3 2 1]))
hold on
plot(ts([p_point(1) p_point(end)])/60/60,[-20 20])
hold off
grid on
title('Surface temp - fitted surface temp')
xlabel('time (hr)')
ylabel('\deltaT (K)')


% figure(117)
% plot(permute(del_nu(1,1,:),[3 2 1]),permute(f(:,p_point,:),[3 1 2]))
% xlabel('wavenumber')
% ylabel('Normalized f')
% grid on

% figure(118)
% plot(permute(del_nu(1,1,:),[3 2 1]),permute(absorption_f(:,p_point,:),[3 1 2]))
% xlabel('wavenumber')
% ylabel('absorption')
% grid on


% figure(777)
% plot(O2_Offline_Backscatter_Channel_Raw_Data(:,1),'--')
% hold on
% plot(O2_Online_Backscatter_Channel_Raw_Data(:,1),'-')
% plot(O2_Offline_Molecular_Detector_Backscatter_Channel_Raw_Data(:,1),'-*')
% plot(O2_Online_Molecular_Detector_Backscatter_Channel_Raw_Data(:,1))
% legend('Off Back','On back','Off mol','On mol')
% hold off


% figure(99587)
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(T_etalon,[3 2 1]))
% grid on
% xlabel('wavenumber')
% title('etalon transmission')

% figure(99588)
% semilogy(permute(del_f(:,p_point,:),[3 2 1]),permute(doppler_O2_ret(:,p_point,:),[3 1 2]))
% hold on
% semilogy(permute(del_f(:,p_point,:),[3 2 1]),permute(T_etalon(:,p_point,:),[3 1 2]))
% hold off
% grid on
% ylim([10^-3 1])
% xlim([-2.5 2.5])
% xlabel('freqency shift (GHz)')
% title('doppler broadending')


% figure(99589)
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(g1(:,p_point,:),[3 1 2]),'-*')
% grid on
% xlabel('wavenumber')
% title('backscatter lineshape')

% figure(99590)
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(eta(:,p_point,:),[3 1 2]),'-*')
% hold on
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:),[3 1 2]),'-*')
% grid on
% xlabel('wavenumber')
% legend('eta','zeta')

% figure(99491)
% plot(W1(:,p_point),rm)
% legend('W1')

figure(99492)
plot(diag(BSR(:,p_point)),rkm,'linewidth',2)
legend('BSR')
title({'Backscatter Ratio';datestr(date_ts_N(p_point(1)))})
grid on
ylabel('Range (km)')

% 
% figure(99493)
% 
%     [~,g] = sgolay(2,oversample);
%     parfor j=1:i_time
%         dBSRSG(:,j) = conv(BSR(:,j), factorial(1)/(-rangeBin)^1 * g(:,2), 'same');
%     end
% plot(diag(dBSRSG(:,p_point)),rkm)
% hold on
% dBSR = (BSR(ind_r_hi,:)-BSR(ind_r_lo,:))/(rangeBin*oversample);
% plot(diag(dBSR(:,p_point)),rkm(ind_r_lo+4))
% 
% dBSRlo = diff(BSR)/rangeBin;
% plot(diag(dBSRlo(:,p_point)),rkm(1:end-1))
% hold off
% legend('2SG','over forward','forward')
% title({'dBSR/dr';datestr(date_ts_N(p_point(1)))})
% grid on
% ylabel('Range (km)')

figure(99494)
plot(diag(alpha_0(:,p_point)-absorption_sonde{t_index}),rkm)

% % % figure(99495)
% % % cla reset
% % % line(diag(dBSRSG(:,p_point)),rkm,'Color','r')
% % % line(diag(dBSRSG(:,p_point)./BSR(:,p_point)),rkm,'Color','b')
% % % %line(diag(dBSRlo(:,p_point)),rkm(1:end-1))
% % % legend('dBSR/dr')
% % % xlim([-.01 .01])
% % % xlabel('dBSR/dr')
% % % ylabel('Range (km)')
% % % 
% % % ax1 = gca; % current axes
% % % ax1.XColor = 'r';
% % % ax1.YColor = 'r';
% % % 
% % % ax1_pos = ax1.Position; % position of first axes
% % % ax2 = axes('Position',ax1_pos,...
% % %     'XAxisLocation','top',...
% % %     'YAxisLocation','right',...
% % %     'Color','none');
% % % 
% % % line(diag(alpha_0(:,p_point)-absorption_sonde{t_index}),rkm,'Parent',ax2,'Color','k')
% % % 
% % % grid on
% % % xlim([-.5 .5]*10^-3)
% % % xlabel('\alpha_0 - \alpha_{sonde}')
% % % ylabel('Range (km)')
% % % legend('\alpha_0 - \alpha_{sonde}')


% figure(99593)
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(Tm0(:,p_point,:),[3 1 2]),'-*')
% grid on
% xlabel('wavenumber')
% title('transmission')

% figure(99594)
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:),[3 1 2]),'-*')
% hold on
% %plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:).*(1-f(:,p_point,:)),[3 1 2]),'-*')
% hold off
% grid on
% xlabel('wavenumber')
% title('zeta1')

% figure(99595)
% %plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:),[3 1 2]),'-*')
% hold on
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:).*(1-f(:,p_point,:)),[3 1 2]),'-*')
% hold off
% grid on
% xlabel('wavenumber')
% title('zeta1s')

% figure(99497)
% plot(zeta_int(:,p_point),rm)
% hold on
% plot(zeta_ls_int(:,p_point),rm)
% hold off
% legend('zeta_int','zeta ls int')
% 
% figure(99498)
% plot(G1(:,p_point),rm)
% legend('G1')

figure(99499)
plot(diag(WV(:,p_point)),rm)
xlabel('WV (molec/m^2')
ylabel('Range (km)')
title({'Water vapor profile';datestr(date_ts_N(p_point(1)))})





% figure(68)
% %p_point = round((8/24)*length(time_Temperature_HSRL));
% plot(Temperature_HSRL(:,p_point),range_Temperature_HSRL)

% figure(465)
% subplot(2,1,1)
% imagesc(thr,rm,cloud_SDm_on)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2on standard deviation')
% subplot(2,1,2)
% imagesc(thr,rm,cloud_SDm_off)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2off standard deviation')

%figure(46)
%plot(rm,o2on_range_corrected(:,p_point))

%figure(47)
%plot(rm,o2on_edge(:,p_point))

% figure(565)
% plot(rm,o2off_SD(:,p_point))
% hold on
% plot(rm,o2on_SD(:,p_point))
% hold off
% legend('on','off')
% title('Standard Deviation filter')
% xlabel('Range(m)')

% figure(56)
% subplot(2,1,1)
% imagesc(thr,rm,o2on_SD)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2on standard deviation')
% subplot(2,1,2)
% imagesc(thr,rm,o2off_SD)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2off standard deviation')

% figure(57)
% subplot(2,1,1)
% imagesc(thr,rm,o2on_edge)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2on edge')
% subplot(2,1,2)
% imagesc(thr,rm,o2off_edge)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2off edge')

% figure(58)
% subplot(2,1,1)
% imagesc(thr,rm,o2on_range_corrected)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2on range corrected')
% subplot(2,1,2)
% imagesc(thr,rm,o2off_range_corrected)
% hold on
% xline(thr(p_point),'r');
% hold off
% set(gca, 'YDir','normal')
% title('o2off range corrected')

% figure(829892)
% plot(o2on_SNR(:,p_point),rm)
% title('o2on_SNR')



figure(884)
plot(diag(T(:,p_point)),rkm,'linewidth',2)
hold on
%%%plot(diag(T(:,p_point)+2),rkm)
%%%plot(diag(T(:,p_point)-2),rkm)
%%%plot(Ts(p_point(1)),0,'+')
%%%plot(Ts_fit(1,p_point(1),end),0,'+')
plot(diag(T_finalm(:,p_point)),rkm,'linewidth',2)

%%%%%%%%%%plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end)),rkm,'--')
%plot(T_final(:,p_point),rm)
exclusion(exclusion==0)=NaN;
%%%plot(diag(T_final_test(:,p_point)).*diag(exclusion(:,p_point,end)),rkm,'*')
plot(T_real(:,t_index),rkm,'.-','linewidth',2)
plot(diag(T_final_test(:,p_point)),rkm)
hold off
grid on
xlabel('Temperature (K)')
ylabel('Range (km)')
ylim([-inf 5])
title({'Temperature profile';datestr(date_ts_N(p_point(1)))})
%legend('Temp guess','Tg+2','Tg-2','Ts','Ts fit','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde','Location','southwest')
legend('Temp guess','Retrieved Temp','Radiosonde Temperature')


figure(984)
plot(diag(P(:,p_point)),rkm)
hold on
%%%plot(diag(T(:,p_point)+2),rkm)
%%%plot(diag(T(:,p_point)-2),rkm)
%%%plot(Ts(p_point(1)),0,'+')
%%%plot(Ts_fit(1,p_point(1),end),0,'+')
plot(diag(Patm_final(:,p_point)),rkm)
%plot(T_final(:,p_point),rkm)
%%%plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end)),rkm,'--')
%plot(T_final(:,p_point),rm)
exclusion(exclusion==0)=NaN;
%%%plot(diag(T_final_test(:,p_point)).*diag(exclusion(:,p_point,end)),rkm,'*')
plot(Patm_real(:,t_index),rkm,'.-')
hold off
grid on
xlabel('Pressure (atm)')
ylabel('Range (km)')
ylim([-inf 5])
title({'Pressure profile';datestr(date_ts_N(p_point(1)))})
%legend('Temp guess','Tg+2','Tg-2','Ts','Ts fit','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde','Location','southwest')
legend('Press guess','Retrieved Press','Radiosonde Press')



% Single time temperature deviation
figure(885)
plot(T_real(:,t_index)-diag(T_finalm(:,p_point)),rkm)
%title(sprintf('Temperature Deviation  %s (UTC)',date_time))
line([0 0],[0 6],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--')
line([-1 -1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.')
line([1 1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
line([2 2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
title({'\Delta T (T_{sonde} - T_{retrieved}) (K)';datestr(date_ts_N(p_point(1)))})
ylabel('Range (km)')
legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','Southwest')
ylim([0 5])
xlim([-10 10])

figure(888)
plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end))-diag(T_finalm(:,p_point)),rkm)
%title(sprintf('Temperature Deviation  %s (UTC)',date_time))
line([0 0],[0 6],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--')
line([-1 -1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.')
line([1 1],[0 4],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
line([2 2],[0 6],'Color',[0.8500, 0.3250, 0.0980])
title({'\Delta T (T_{fit} - T_{retrieved}) (K)';datestr(date_ts_N(p_point(1)))})
ylabel('Range (km)')
legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','Southwest')
ylim([0 5])
xlim([-10 10])

% figure(9985)
% semilogy(thr,bg_o2off)
% hold on
% semilogy(thr,bg_o2on)
% semilogy(thr,bg_o2off_mol)
% semilogy(thr,bg_o2on_mol)
% title('Background')
% xlabel('Time UTC hr')
% ylabel('photons')
% xline(thr(p_point));
% hold off
% legend('Offline','Online','Offline molecular','Online molecular')
% grid on


figure(9985)
bg_o2offm = bg_o2off.*mask_nodata(1,:);
bg_o2offm(bg_o2offm <= 0) = NaN;
[bg_o2offm,TF]=rmmissing(bg_o2offm,2);
bg_o2onm = bg_o2on.*mask_nodata(1,:);
bg_o2onm(bg_o2onm <= 0) = NaN;
bg_o2onm=rmmissing(bg_o2onm,2);
bg_o2off_molm = bg_o2off_mol.*mask_nodata(1,:);
bg_o2off_molm(bg_o2off_molm <= 0) = NaN;
bg_o2off_molm=rmmissing(bg_o2off_molm,2);
bg_o2on_molm = bg_o2on_mol.*mask_nodata(1,:);
bg_o2on_molm(bg_o2on_molm <= 0) = NaN;
bg_o2on_molm=rmmissing(bg_o2on_molm,2);

semilogy(datenum(date_ts_N(~TF)),bg_o2offm)
hold on
semilogy(datenum(date_ts_N(~TF)),bg_o2onm)
semilogy(datenum(date_ts_N(~TF)),bg_o2off_molm)
semilogy(datenum(date_ts_N(~TF)),bg_o2on_molm)
legend('Offline','Online','Offline molecular','Online molecular')
xlim([datenum(date_ts_N(1)) datenum(date_ts_N(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
xtickangle(tickAngle)
title('Background')
xlabel('Time UTC')
ylabel('photons')
hold off
grid on

% figure(66666)
% surf(ts,rm,T_finalm)
% title('retrieved temperature')

% Full retrieved temperature profile with SNR mask
% figure(885)
% r_max_plot = 5; % Max range to plot [km]
% [~,ind_km_max] = min(abs(rkm-r_max_plot));
% 
% imAlpha=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
% 
% colorLimits = [min(T_finalm(1:ind_km_max,:),[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')];
% imagesc(thr,rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% set(gca, 'YDir','normal')
% set(gca,'color',0*[1 1 1]);%color background black
% colorbar
% colormap(parula)
% hold on
% 
% imAlpha2=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha2(cloud_SDm(1:ind_km_max,:).*cloud_SDm_off(1:ind_km_max,:) == 1)=0;
% mapCloud = imagesc(thr,rkm(1:ind_km_max),cloud_SDm(1:ind_km_max,:).*cloud_SDm_off(1:ind_km_max,:),'AlphaData',imAlpha2)
% colormap(gca,white)
% 
% xline(thr(p_point),'r')
% hold off
% %title(sprintf('Temperature %s to %s (UTC)',ts_datetime(1),ts_datetime(end)))
% xlabel('Time (hours)')
% ylabel('Range (km)')
% title(colorbar,'Temperature (K)')


%    % New figure
%     hf=figure(887);
%     r_max_plot = 5; % Max range to plot [km]
% [~,ind_km_max] = min(abs(rkm-r_max_plot));
% 
% imAlpha=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
% 
% colorLimits = [min(T_finalm(1:ind_km_max,:),[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')];
%     %Background image
%     h1 = axes;colormap(h1,'parula');
%     p1=imagesc(thr,rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits);
%     set(h1,'ydir','normal');
%     set(gca,'color',[1 1 1]);%color background black
%     colorbar(h1)
%     %Foreground image
%     hold on
%     
%     imAlpha2=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha2(cloud_SDm(1:ind_km_max,:).*cloud_SDm_off(1:ind_km_max,:) == 1)=0;
%     h2=axes;
%     p2=imagesc(thr,rkm(1:ind_km_max),cloud_SDm(1:ind_km_max,:).*cloud_SDm_off(1:ind_km_max,:),'AlphaData',imAlpha2);
%     set(h2,'color','black','visible','off')
%     colormap(h2,'gray');
%     set(h2,'ydir','normal');
% 
%     hold off
%     colorbar(h2,'hide')
%     linkaxes([h1 h2])


% figure(886)
% % if length(span_days)==1
% %     title_temp = append('Temperature profile for ' , datestr(span_days));
% % else
% %     title_temp = 'Temperature profile for ';
% %     for i = 1:length(span_days)
% %         title_temp = append(title_temp , datestr(span_days(i)));
% %     end
% % end
% r_max_plot = 5; % Max range to plot [km]
% [~,ind_km_max] = min(abs(rkm-r_max_plot));
% imAlpha=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
% imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
% %imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
% colorLimits = [min(T_finalm(1:ind_km_max,:)-1,[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')];
% imagesc(thr,rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% hold on
% xline(thr(p_point),'color','b');
% colors = colormap;
% colors(1,:) = [0 0 0];%set lowest color black
% colormap(colors);
% set(gca, 'YDir','normal')
% set(gca,'color',[1 1 1]);%color background white
% colorbar
% ylim([-inf 4])
% %title(sprintf('Temperature %s to %s (UTC)',ts_datetime(1),ts_datetime(end)))
% xlabel('Time (UTC hours)')
% ylabel('Range (km)')
% title(colorbar,'Temperature (K)')
% hold off
% title('Temperature retrieval on 04/19/2019')
% %title(file)


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
imAlpha=ones(size(T_finalm(1:ind_km_max,:)));
imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(T_finalm(1:ind_km_max,:)-1,[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(1,i))
        %disp('ran')
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end

xline(datenum(date_ts(p_point(1))))
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('Temperature %s to %s (UTC)',date_ts(1),date_ts(end)))
%title(sprintf('Temperature retrieval\n%s(UTC)',span_days(1)))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'Temperature (K)')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off
%title('Temperature retrieval')
%title(file)


figure(8886)
subplot(5,1,1)
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(alpha_0m(1:ind_km_max,:)));
imAlpha(isnan(alpha_0m(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_0m(1:ind_km_max,:),[],'all'), max(alpha_0m(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),alpha_0m(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(1,i))
        %disp('ran')
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end
xline(datenum(date_ts(p_point(1))))

colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha 0 %s to %s (UTC)',date_ts(1),date_ts(end)))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off


subplot(5,1,2)
alpha_1m = alpha_1 .* SNRm .* cloud_SDm_above;
%%alpha_1m(alpha_1m <= 0) = NaN;          % Replace mask with NaNs
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(alpha_1m(1:ind_km_max,:)));
imAlpha(isnan(alpha_1m(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_1m(1:ind_km_max,:),[],'all'), max(alpha_1m(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),alpha_1m(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(1,i))
        %disp('ran')
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end
xline(datenum(date_ts(p_point(1))))
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha 1'))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off


subplot(5,1,3)
alpha_2m = alpha_2 .* SNRm .* cloud_SDm_above;
%%alpha_2m(alpha_2m <= 0) = NaN;          % Replace mask with NaNs
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(alpha_2m(1:ind_km_max,:)));
imAlpha(isnan(alpha_2m(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_2m(1:ind_km_max,:),[],'all'), max(alpha_2m(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),alpha_2m(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(1,i))
        %disp('ran')
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end
xline(datenum(date_ts(p_point(1))))
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha 2'))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off


subplot(5,1,4)
alpha_total_rawm = alpha_total_raw .* SNRm .* cloud_SDm_above;
%%alpha_2m(alpha_2m <= 0) = NaN;          % Replace mask with NaNs
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(alpha_total_rawm(1:ind_km_max,:)));
imAlpha(isnan(alpha_total_rawm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_total_rawm(1:ind_km_max,:),[],'all'), max(alpha_total_rawm(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),alpha_total_rawm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(1,i))
        %disp('ran')
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end
xline(datenum(date_ts(p_point(1))))
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha raw'))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off


subplot(5,1,5)
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(alpha_totalm(1:ind_km_max,:)));
imAlpha(isnan(alpha_totalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
%imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_totalm(1:ind_km_max,:),[],'all'), max(alpha_totalm(1:ind_km_max,:),[],'all')];
imagesc(datenum(date_ts),rkm(1:ind_km_max),alpha_totalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on

%Plotting sode locations
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(1,i))
        %disp('ran')
        %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
        plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
    end
end
xline(datenum(date_ts(p_point(1))))
    
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha total'))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
%xticklabels(datestr(date_ts))
xlim([datenum(date_ts(1)) datenum(date_ts(end))]) 
xticks(datenum(date_ts(1)): tickSpacing :datenum(date_ts(end)))
xtickangle(tickAngle)
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
hold off



% figure(888)
% % if length(span_days)==1
% %     title_temp = append('Temperature profile for ' , datestr(span_days));
% % else
% %     title_temp = 'Temperature profile for ';
% %     for i = 1:length(span_days)
% %         title_temp = append(title_temp , datestr(span_days(i)));
% %     end
% % end
% r_max_plot = 5; % Max range to plot [km]
% [~,ind_km_max] = min(abs(rkm-r_max_plot));
% imAlpha=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
% imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
% colorLimits = [min(T_finalm(1:ind_km_max,:)-1,[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')];
% contourf(thr,rkm(1:ind_km_max),T_finalm(1:ind_km_max,:));%,'AlphaData',imAlpha,colorLimits)
% colors = colormap;
% colors(1,:) = [0 0 0];%set lowest color black
% colormap(colors);
% set(gca, 'YDir','normal')
% set(gca,'color',[1 1 1]);%color background black
% colorbar
% hold on
% imAlpha2=ones(size(T_finalm(1:ind_km_max,:)));
% imAlpha2(cloud_SDm(1:ind_km_max,:) == 1)=0;
% imagesc(thr,rkm(1:ind_km_max),cloud_SDm(1:ind_km_max,:),'AlphaData',imAlpha2)
% xline(thr(p_point),'color','b');
% hold off
% %title(sprintf('Temperature %s to %s (UTC)',ts_datetime(1),ts_datetime(end)))
% xlabel('Time (hours)')
% ylabel('Range (km)')
% title(colorbar,'Temperature (K)')
% % title(title_temp)

% figure(2)
% surf(ts,rm,o2on)
% title('Averaged return photons online')
% xlabel('Time s')
% ylabel('Range m')
% zlabel('Photons')
% 
% figure(3)
% surf(ts,rm,o2off)
% title('Averaged return photons offline')
% xlabel('Time s')
% ylabel('Range m')
% zlabel('Photons')



% figure(5)
% surf(thr,rm/1000,alpha_0m)
% xlabel('Relative time [hours]')
% ylabel('range [km]')
% zlabel('Zeroth order absorption[m-1]')

figure(4)
plot(diag(alpha_0(:,p_point)),rkm,'linewidth',2)
hold on
%%plot(diag(absorption(:,p_point)),rkm)
%plot(permute(absorption_f(:,p_point,online_index(1)-2:online_index(1)+2),[1 3 2]),rkm,'-*')
%%%%plot(permute(absorption_f(:,p_point,online_index(1)),[1 3 2]),rkm,'-*')
%plot(rm,absorption_new(:,online_index_new)-absorption_new(:,offline_index_new))
plot(diag(alpha_totalm(:,p_point)),rkm,'linewidth',2)
%%plot(diag(alpha_0m(:,p_point))+diag(alpha_1(:,p_point)),rkm)
%%plot(diag(alpha_total_raw(:,p_point)),rkm)
plot(absorption_sonde{t_index},rkm,'-.','linewidth',2)
plot(diag(alpha_0(:,p_point)+alpha_1(:,p_point)),rkm)
plot(diag(alpha_0(:,p_point)+alpha_1(:,p_point)+alpha_2(:,p_point)),rkm)

%plot(diag(alpha_02(:,p_point)),rkm-.150/2,'.-')
hold off
grid on
ylim([0 5])
xlim([0 5]*10^-4)
title({'O2 absorption profile';datestr(date_ts_N(p_point(1)))})
ylabel('Range (km)')
xlabel('Absorption m^{-1}')
%legend('\alpha_0','\alpha_{total}','\alpha_{sonde}','\alpha_0 + \alpha_1','\alpha_0 + \alpha_1 + \alpha_2')
%legend('\alpha_0','\alpha_{total}','\alpha_{sonde}','\alpha_0 + \alpha_1 + \alpha_2')
%legend('\alpha_0','absorption','\alpha_{total}','\alpha_{sonde}','\alpha_0 + \alpha_1 + \alpha_2')
%legend('Measured zeroth order absorption','Theoretical absorption from Tg','Pertabative absorption','Alpha 0+1','Alpha 0+1+2','Sonde')%,'theoretical absorption online-offine')
legend('Measured zeroth order absorption','Corrected absorption','Radiosonde absorption')%,'theoretical absorption online-offine')


% figure(6)
% plot(rm,alpha_0m(:,round(end/2)))
% title('Calculated absorption coefficient with SNR mask')
% xlabel('Range m')
% ylabel('Absorption m^-1')

% figure(7)
% plot(rm,o2on(:,p_point))
% hold on
% plot(rm,o2off(:,p_point))
% title('Averaged photon returns')
% legend('Online','Offline')
% hold off
% grid on
% xlabel('Range [m]')
% ylabel('Photons')

figure(8)
plot(rm,o2on(:,p_point(1)).*rm.^2)
hold on
plot(rm,o2off(:,p_point(1)).*rm.^2)
%title('Averaged photon returns multiplied by r^2')
title({'Averaged counts returns multiplied by r^2';datestr(date_ts_N(p_point(1)))})
legend('Online','Offline')
xlabel('Range m')
ylabel('Counts * r^2')
grid on
hold off


figure(7)
i_range = length(rm);

plot(o2on(:,p_point(1)),rm)
hold on
plot(o2off(:,p_point(1)),rm)
plot(o2on_mol(:,p_point(1)),rm)
plot(o2off_mol(:,p_point(1)),rm)

plot(o2on_noise(1:i_range,p_point(1)),rm,'--')
plot(o2off_noise(1:i_range,p_point(1)),rm,'--')
plot(o2on_noise_mol(1:i_range,p_point(1)),rm,'--')
plot(o2off_noise_mol(1:i_range,p_point(1)),rm,'--')
title({'Averaged counts';datestr(date_ts_N(p_point(1)))})
legend('Online','Offline','Online molecular','Offline molecular')
hold off
grid on
ylabel('Range [m]')
xlabel('Counts')
xlim([-inf 300])

% figure(77)
% Y = fft(o2on_noise(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2on = P2(1:i_range/2+1);
% 
% Y = fft(o2off_noise(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2off = P2(1:i_range/2+1);
% 
% Y = fft(o2on_noise_mol(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2on_mol = P2(1:i_range/2+1);
% 
% Y = fft(o2off_noise_mol(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2off_mol = P2(1:i_range/2+1);
% 
% Y = fft(o2on(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2on2 = P2(1:i_range/2+1);
% 
% Y = fft(o2off(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2off2 = P2(1:i_range/2+1);
% 
% Y = fft(o2on_mol(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2on_mol2 = P2(1:i_range/2+1);
% 
% Y = fft(o2off_mol(1:i_range,p_point(1)));
% P2 = abs(Y/i_range);
% P1o2off_mol2 = P2(1:i_range/2+1);
% 
% krange = [zeros(1,62) (1/8)*ones(1,8) zeros(1,63)];
% krange = [zeros(1,(i_range-oversample)/2) (1/oversample)*ones(1,oversample) zeros(1,(i_range-oversample)/2)];
% %krange = [zeros(1,62) ones(1,8) zeros(1,63)];
% Y = fft(krange');
% P2 = abs(Y/i_range);
% P1k = P2(1:i_range/2+1);
% 
% ksg = sgolay(0,oversample);
% ksg = ksg((oversample+1)/2,:);
% ksg = [zeros(1,62) ksg zeros(1,62)];
% Y = fft(ksg');
% P2 = abs(Y/i_range);
% P1ksg = P2(1:i_range/2+1);
% 
% f = (1/rangeBin)*(0:(i_range/2))/i_range;
% 
% semilogy(f,P1o2on)
% hold on
% semilogy(f,P1o2off)
% semilogy(f,P1o2on_mol)
% semilogy(f,P1o2off_mol)
% 
% semilogy(f,P1o2on2,'.-')
% semilogy(f,P1o2off2,'.-')
% semilogy(f,P1o2on_mol2,'.-')
% semilogy(f,P1o2off_mol2,'.-')
% 
% semilogy(f,P1k,'--')
% semilogy(f,P1ksg,'--')
% semilogy(f,P1k.*P1o2on,'--')
% hold off
% grid on
% legend('on','off','on mol','off mol')
% title('Single sided spectrum of o2on')
% xlabel('f (1/m)')
% ylabel('|P1(f)|')


figure(777)
Y = fft2(o2on);
imagesc(log(abs(fftshift(Y))));

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

if ~isempty(data_col_real)

figure(6773)
%temp diff sonde
hold on
index=1;
for i = 1:size(data_col_real,2)
    if ~isnan(data_col_real(i))
        plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*')
        fullTFinal(:,index) = T_finalm(:,data_col_real(i));
        T_real_final(:,index) = T_real(:,i);
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

figure(6774)
for i = 1:length(270:305)
    for j = 1:length(270:305)
        
        tempProb(i,j) = nnz((fullTFinal<(270+i-1)+0.5)&(fullTFinal>=(270+i-1)-0.5)& (T_real_final<(270+j-1)+0.5)&(T_real_final>=(270+j-1)-0.5));
    end
end
imagesc(270:305,270:305,tempProb')
colorbar
set(gca, 'YDir','normal')
xlabel('Temp retrieval (K)')
ylabel('Sonde (K)')



figure(6775)
cloud_index = [6 11 12 13 14];
cloud_index = [1 3 4 5 7 8 9];
cloud_index = [6 7 8];
cloud_index = [4];
%temp diff sonde
hold on
index=1;
for i = cloud_index
    if ~isnan(data_col_real(i))
        plot(T_finalm(:,data_col_real(i)),T_real(:,i),'*')
        fullTFinal(:,index) = T_finalm(:,data_col_real(i));
        T_real_final(:,index) = T_real(:,i);
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

end

BSRm = BSR .* SNRm .* cloud_SDm .* mask_nodata;
BSRm(BSRm <= 0) = NaN;


figure(7473)
imagesc(thr,rm,BSRm)
hold on
xline(thr(p_point(1)),'r');
hold off
colorbar
colormap(flipud(hot))
caxis([0 5])
set(gca, 'YDir','normal')
title('BSR')


figure(589)
subplot(2,1,1)
imagesc(thr,rm/1000,o2on_noise.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
%caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')
title('o2 on unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,1,2)
imagesc(thr,rm,o2off_noise.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
%set(gca,'ColorScale','log')
colorbar
title('o2off unaveraged')



