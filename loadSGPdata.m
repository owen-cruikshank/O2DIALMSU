function [Range,Time,Counts,Sonde,Model,Spectrum,BSR,Data] = loadSGPdata(span_days,Options,Constant)

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
% % addpath('C:\Users\oencr\OneDrive - Montana State University\Research\O2 DIAL\Data\ARM DOE SGP Data')
% % addpath('C:\Users\d98b386\OneDrive - Montana State University\Research\O2 DIAL\Data\ARM DOE SGP Data')

%%
analysispath = pwd;

for i = 1:length(span_days)
 date_i = span_days(i);
    d = num2str(yyyymmdd(date_i));
    filenameMPD = insertAfter('wv_dial05..Python.nc',10,d(3:end));
    %filenameMPD = insertAfter('mpd05..Python.nc',6,d(3:end));
    
    cd ..
    cd('Data\ARM DOE SGP Data\')
   filenameMPD = fullfile(pwd,filenameMPD);
    
    filenameSGP = dir(['sgpsondewnpnC1.b1.',d,'.*']);
%     for ii=1:length(filenameSGP)
%     filenameSGP{ii} = fullfile(pwd,filenameSGP{ii});
%     end
    
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
        [T_sgp_raw{j,i},P_sgp_raw{j,i},rm_sgp{j,i},t_single_sgp{j,i}] = getSGPdata(fullfile(filenameSGP(j).folder,filenameSGP(j).name));
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

file = filenameMPD;
%file = ['ARM DOE SGP Data\mpd05.' d(3:end) '.Python.nc'];
%file = ['mpd05.' d(3:end) '.Python.nc'];

%read in relevant data
time_O2_Online_Backscatter_Channel_Raw_Data{i}                             = ncread(file,'time_O2_Online_Backscatter_Channel_Raw_Data');                           %[s](1440x1){single} seconds since midnight
time_O2_Online_Backscatter_Channel_Raw_Data{i} = time_O2_Online_Backscatter_Channel_Raw_Data{i} + (i-1)*86400;
range_O2_Online_Backscatter_Channel_Raw_Data{i}                            = ncread(file,'range_O2_Online_Backscatter_Channel_Raw_Data');                          %[m](560x1){single}
O2_Online_Backscatter_Channel_Raw_Data{i}                                  = ncread(file,'O2_Online_Backscatter_Channel_Raw_Data');                                %[photons](560x14400){single}Unpolarization Online O2-DIAL Backscatter Returns on the Combined Detector'
%O2_Online_Backscatter_Channel_Raw_Data_variance                         = ncread(file,'O2_Online_Backscatter_Channel_Raw_Data_variance');
O2_Online_Backscatter_Channel_Raw_Data_ProfileCount{i}                     = ncread(file,'O2_Online_Backscatter_Channel_Raw_Data_ProfileCount');                   %[count](1440x1){int}Number of raw profiles integrated into this profile
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
O2_Online_Backscatter_Channel_Raw_Data_ProfileCount = double(cat(1,O2_Online_Backscatter_Channel_Raw_Data_ProfileCount{:}));
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

Counts.NBins = O2_Online_Backscatter_Channel_Raw_Data_ProfileCount';

%convert time from single to double
ts_raw_o2_on = double(time_O2_Offline_Backscatter_Channel_Raw_Data);   %[s] Raw time vector from data files 
ts_raw_o2_off = double(time_O2_Online_Backscatter_Channel_Raw_Data);   %[s] Raw time vector from data files 

o2on_raw = double(O2_Online_Backscatter_Channel_Raw_Data);      
o2off_raw = double(O2_Offline_Backscatter_Channel_Raw_Data);   


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
date_ts_N = datetime(y,m,d,0,0,Time.ts);
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
Range.NBins = double(560);                                %[none] Convert data to double
Range.rangeBin = (Constant.c * Range.nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
%rangeBin = rangeBin*2;
Range.rm_raw_o2 = Range.rangeBin:Range.rangeBin:Range.NBins(1)*Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector
Range.r_max = 5000;                                       %[m] Max range 
Range.rm = Range.rm_raw_o2(Range.rm_raw_o2<=Range.r_max);                   %[m] Shorten range vector
Range.rkm = Range.rm./1000;                                     %[km] Range vector
Range.i_range = length(Range.rm);                               %[none] Size of range vector

%%
%%
% SGP Radiosonde Data
% Can't interpolate with array with repeating values, so use Kevin's sonde preparation function

increment = 1;
for i = 1:numel(T_sgp_raw)
    if ~isempty(T_sgp_raw{i}) & ~isnan(T_sgp_raw{i})
        % Subtract first range value (site elevation) from whole vector
        rm_sgp{i} = rm_sgp{i} - rm_sgp{i}(1);
        % Collect radiosonde surface measurements
        T_sgp_surf(increment) = T_sgp_raw{i}(1);
        P_sgp_surf(increment) = P_sgp_raw{i}(1);
        % Convert datetimes from cells to vector
        datetime_sgp(increment) = datetime_sgp_cell{i};
        % Custom interpolation function
        [T_sgp_int{increment},P_sgp_int{increment},rm_sgp_int{increment}] = interp_sonde(T_sgp_raw{i},P_sgp_raw{i},rm_sgp{i},Range.rangeBin);
        T_sgp(:,increment) = T_sgp_int{increment}(1:length(Range.rm));
        P_sgp(:,increment) = P_sgp_int{increment}(1:length(Range.rm));
        
        %interp sonde time
        [rm_sgp{i},IA,IC] = unique(rm_sgp{i});
        
        sonde_time(1:length(Range.rm),increment) = interp1(rm_sgp{i},t_single_sgp{i}(IA),Range.rm)';
        %sonde_time = sonde_time';
        
                 %Find index of sonde in time vector
         for j = 1:length((Range.rm))
            [~, Sonde.sonde_ind(j,increment)]=min(abs(sonde_time(j,increment)-date_ts_N));
         end
         
         increment=increment+1;
    else
%         T_sgp(1:length(Range.rm),i) = nan(length(Range.rm),1);
%         P_sgp(1:length(Range.rm),i) = nan(length(Range.rm),1);
%         datetime_sgp(1:length(Range.rm),i)=nan(length(Range.rm),1);
%         sonde_time(1:length(Range.rm),i) = nan(length(Range.rm),1);
%         Sonde.sonde_ind(1:length(Range.rm),i) = nan(length(Range.rm),1);
        
    end
end



Sonde.T_sonde = T_sgp;
Sonde.P_sonde = P_sgp;
datetime_real = datetime_sgp;
i_time_real = length(datetime_real);

% Testing to see if starting at 60 m will help.... ************************************************************************
%%%T_surf = T_60m;
%%%P_surf = P_60m;

% Find time lapsed since beginning datetime
time_lapsed = datetime_real - span_days(1);
[hrs_real,mins_real,sec_real] = hms(time_lapsed);

% Indices where "real" data is available
%data_col_real = hrs_real*60 + mins_real + 1;

%%%%%[~,data_col_real] = min(abs(ts' - (hrs_real*60*60+mins_real*60)));

%%
% Count dead time correction
nsPerBin = 250e-9;
summedBins = 14000/2;
%summedBins = 1000/2;
%summedBins = 560;
% Deadtime for apd SPCM-780-12-FC
deadTime = 22e-9;
o2on_countRate = o2on_raw /(nsPerBin * summedBins);
o2on_correctionFactor = 1./(1 - deadTime*o2on_countRate);
%o2on_correctionFactor=1;
o2on_raw_corrected = o2on_correctionFactor .* o2on_raw;
o2off_countRate = o2off_raw /(nsPerBin * summedBins);
o2off_correctionFactor = 1./(1 - deadTime*o2off_countRate);
%o2off_correctionFactor=1;
o2off_raw_corrected = o2off_correctionFactor .* o2off_raw;

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
% o2on_intp = interp2(Time.ts,rm_raw,o2on_intp,Time.ts,Range.rm_raw_o2,'nearest',nan);
% o2off_intp = interp2(Time.ts,rm_raw,o2off_intp,Time.ts,Range.rm_raw_o2,'nearest',nan);
% o2on_intp_mol = interp2(Time.ts,rm_raw,o2on_intp_mol,Time.ts,Range.rm_raw_o2,'nearest',nan);
% o2off_intp_mol = interp2(Time.ts,rm_raw,o2off_intp_mol,Time.ts,Range.rm_raw_o2,'nearest',nan);

%%
% ======================
% === O2 Backscatter ===
% ======================
% --- O2 Background Subtraction ---
% Online
Counts.bg_o2on = mean(o2on_intp(end-150:end,:),1,'omitnan');% Take mean of last data points
Counts.o2on_bgsub = o2on_intp - Counts.bg_o2on;       % Background subtracted
Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
Counts.o2on_noise = Counts.o2on_bgsub(1:Range.i_range,:);   % Shorten to length of range vector

% Offline
Counts.bg_o2off = mean(o2off_intp(end-150:end,:),1,'omitnan');% Take mean of last data points
Counts.o2off_bgsub = o2off_intp - Counts.bg_o2off;      % Background subtracted
Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
Counts.o2off_noise = Counts.o2off_bgsub(1:Range.i_range,:);   % Shorten to length of range vector

% Online
Counts.bg_o2on_mol = mean(o2on_intp_mol(end-150:end,:),1,'omitnan');% Take mean of last data points
Counts.o2on_bgsub_mol = o2on_intp_mol - Counts.bg_o2on_mol;       % Background subtracted
Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
Counts.o2on_noise_mol = Counts.o2on_bgsub_mol(1:Range.i_range,:);   % Shorten to length of range vector

% Offline
Counts.bg_o2off_mol = mean(o2off_intp_mol(end-150:end,:),1,'omitnan');% Take mean of last data points
Counts.o2off_bgsub_mol = o2off_intp_mol - Counts.bg_o2off_mol;      % Background subtracted
Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero
Counts.o2off_noise_mol = Counts.o2off_bgsub_mol(1:Range.i_range,:);   % Shorten to length of range vector



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
% load('OverlapONOFF.mat','OverlapONOFF')
% Counts.o2off = Counts.o2off.*OverlapONOFF;

% load('overlapCorrSGP061821.mat','overlapCorr')
% %Counts.o2off = Counts.o2off.*overlapCorr;
% [~,constAlt]=min(abs(2000-Range.rm));
% Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*overlapCorr(1:constAlt,:);


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
    BSR(:,ii) = interp1(double(range_Backscatter_Ratio),BSRintSum(:,ii),Range.rm);
end
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
if ~isempty(sonde_time)
    for i=1:numel(sonde_time(1,:))
            if isdatetime(sonde_time(1,i))
                %absorption_sonde{i} = diag(absorption_O2_770_model(T_real(:,i),Patm_real(:,i),nu_online(data_col_real(:,i)),WV(:,data_col_real(:,i)))); %[m-1] Function to calculate theoretical absorption
                Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(1),Model.WV(:,Sonde.sonde_ind(:,i)))); %[m-1] Function to calculate theoretical absorption
                Sonde.trasmission_sonde{i} = exp(-cumtrapz(Range.rm,Sonde.absorption_sonde{i})); %O2 transmission
                
            else
                Sonde.absorption_sonde{i} = nan(Range.i_range,1);
                Sonde.trasmission_sonde{i} = nan(Range.i_range,1);
                Sonde.T_sonde = nan(Range.i_range,1);
                Sonde.P_sonde = nan(Range.i_range,1);
            end
    end
else
      Sonde.absorption_sonde{1} = nan(i_range,1);
      Sonde.T_sonde = nan(i_range,1);
      Sonde.P_sonde = nan(i_range,1);
      
%       absorption_sonde = {};
%       T_sonde = [];
%       P_sonde = [];
end
%%
Data.Thermocouple.InsideCell.TimeStamp = [];
Data.Thermocouple.InsideCell.Temperature = [];
Data.Thermocouple.OutsideCell.Temperature = [];
Data.Thermocouple.TSOA.Temperature = [];
Data.Thermocouple.RoomTemp.Temperature = [];
Data.Laser.O2Online.TimeStamp = [];
Data.Laser.O2Online.WaveDiff = [];
Data.Laser.O2Offline.WaveDiff = [];
Data.Etalon.O2Etalon.TimeStamp = [];
Data.Laser.O2Online.TemperatureActual = [];
Data.Laser.O2Online.TemperatureDesired = [];
Data.Laser.O2Offline.TemperatureActual = [];
Data.Laser.O2Offline.TemperatureDesired = [];
Data.Laser.TWSOA.TemperatureActual = [];
Data.Laser.TWSOA.TemperatureDesired = [];
Data.Etalon.O2Etalon.TemperatureActual = [];
Data.UPS.all.TimeStamp = [];
Data.UPS.all.BatteryTimeLeft = [];
Data.UPS.all.HoursOnBattery = [];
Data.UPS.all.UPSTemperature = [];

Data.MCS.Channel0.Data = nan(length(Range.rm_raw_o2),length(Time.ts));
Data.MCS.Channel2.Data = nan(length(Range.rm_raw_o2),length(Time.ts));
Data.MCS.Channel8.Data = nan(length(Range.rm_raw_o2),length(Time.ts));
Data.MCS.Channel10.Data = nan(length(Range.rm_raw_o2),length(Time.ts));

Data.Laser.TWSOA.TemperatureActual = [];
Data.Laser.TWSOA.TemperatureDesired = [];
Data.Laser.TWSOA.TimeStamp = [];
Data.Laser.O2Online.SeedPower =[];
Data.Laser.O2Offline.SeedPower = [];
Data.Laser.O2Offline.TimeStamp = [];

%Reset Counts
Counts.NBins = o2onNBins;




