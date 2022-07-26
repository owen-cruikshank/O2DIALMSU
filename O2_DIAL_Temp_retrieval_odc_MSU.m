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
m_air = 4.792E-26;          %[kg] Mass of air
q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)
%%

disp('Reading in files')
 
%Read data from MSU DIAL
date_start = datetime(2021,3,28,'TimeZone','UTC');%yyyy,mm,dd
date_start = datetime(2021,4,3,'TimeZone','UTC');%yyyy,mm,dd
date_end = date_start;
%date_end   = datetime(2020,10,8,'TimeZone','UTC');%yyyy,mm,dd
span_days = date_start:date_end;


cd ../
path = [pwd '\Data\'];
cd ../../../
homepath = pwd;
sondepath = [pwd '\Box\Radiosondes\Data\All Data\'];


%sondepath = 'C:\Users\d98b386\Box\Radiosondes\Data\All Data\';
%sondepath = 'C:\Users\oencr\Box\Radiosondes\Data\All Data\';

%cd( '.\OneDrive - Montana State University - Bozeman\Research\O2 DIAL\analysis')
cd( '.\OneDrive - Montana State University\Research\O2 DIAL\analysis')

% path = 'C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\o2DIAL_data';
% %path = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\o2DIAL_data';
% %path = 'C:\Users\Owen C\OneDrive - Montana State University - Bozeman\research s19\o2DIAL_data';
% 
% sondepath = 'C:\Users\d98b386\Box\Radiosondes\Data\All Data\';
% %sondepath = 'D:\Owen\Box\Radiosondes\Data\All Data\';
% %sondepath = 'C:\Users\Owen C\Box\Radiosondes\Data\All Data\';

%%
%load data

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

date_ts = date_start + seconds(ts);


%t_pulse = 1;                    %[microseconds] Pulse duration

% =====================
% === Range Vectors ===
% =====================
% Create range vector from length of data bins
nsPerBin = double(250);                             %[ns] Nanoseconds per bin
NBins = double(560);                                %[none] Number of range bins
rangeBin = (c * nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
rangeMin = -150; %[m]
rangeMin = -300; %[m]
rangeMin = -rangeBin; %[m]
rangeMin = -(c * (1*10^-6))/2;
rangeMin = 0;
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
disp('Loading weather Station Data')

%[weather_Temperature_interp, weather_absPressure_interp, weather_WV_interp] = ORSLweatherv2(span_days,ts,path);
[weather_Temperature_interp, weather_absPressure_interp, weather_WV_interp] = ORSLweatherv3(span_days,ts,path);

% [weather_Temperature_interp, weather_absPressure_interp, weather_WV_interp] = wunderWeather(span_days,ts,path);
% weather_Temperature_interp = weather_Temperature_interp+10;

%Use if there is no data from weather station
% weather_Temperature_interp = ones(1,length(ts))*(17);
% weather_absPressure_interp = ones(1,length(ts))*1013.25;
% weather_absPressure_interp = ones(1,length(ts))*850;
 weather_WV_interp = zeros(1,length(ts));
 
 %weather_Temperature_interp = 272.310000000000*ones(1,length(ts)) -273.15;
 %weather_absPressure_interp = 0.836910930175179*ones(1,length(ts)).* 1013.25;

%old weather
%%%%%%%%[weather_Temperature_interp, weather_absPressure_interp, weather_WV_interp] = ORSLweather(span_days,ts,path);

%=============================
%== Cobleigh sonde =========
%===========================
disp('Loading Sonde data')

[sonde_datetime,sondeStruc] =  COBradiosonde(sondepath,span_days);


for i = 1:numel(sonde_datetime)
    if isdatetime(sonde_datetime(i))%~isnan(sonde_datetime)

        
        %sondeStruc.T = sondeStruc.T + (Ts(sonde_ind)-sondeStruc.T(1));
        % Subtract first range value (site elevation) from whole vector
        rm_sgp{i} = sondeStruc(i).Height - sondeStruc(i).Height(1);
        %convert to same units
        sondeStruc(i).P = sondeStruc(i).P./1013.25;%atm
        % Collect radiosonde surface measurements
        T_sgp_surf(i) = sondeStruc(i).T(1);
        P_sgp_surf(i) = sondeStruc(i).P(1);
        % Convert datetimes from cells to vector
        %datetime_sgp(i) = datetime_sgp_cell{i};
        % Custom interpolation function
        %[T_sonde_int{i},P_sonde_int{i},WV_sonde_int{i},rm_sonde_int{i}] = interp_sonde(sondeStruc(i).T,sondeStruc(i).P,sondeStruc(i).WV,rm_sgp{i},rangeBin);
        [T_sonde_int{i},P_sonde_int{i},WV_sonde_int{i},rm_sonde_int{i}] = interp_sonde2(sondeStruc(i).T,sondeStruc(i).P,sondeStruc(i).WV,rm_sgp{i},rangeBin);
        
        if length(T_sonde_int{i})<length(rm)
            disp('ran')
            T_sonde(1:length(T_sonde_int{i}),i) = T_sonde_int{i};
            T_sonde(length(T_sonde_int{i})+1:length(rm),i)=nan(length(rm)-length(T_sonde_int{i}),1);
            P_sonde(1:length(T_sonde_int{i}),i) = P_sonde_int{i};
            P_sonde(length(P_sonde_int{i})+1:length(rm),i)=nan(length(rm)-length(P_sonde_int{i}),1);
            WV_sonde(1:length(T_sonde_int{i}),i) = WV_sonde_int{i};
            WV_sonde(length(WV_sonde_int{i})+1:length(rm),i)=nan(length(rm)-length(WV_sonde_int{i}),1);
        else
            T_sonde(:,i) = T_sonde_int{i}(1:length(rm));
            P_sonde(:,i) = P_sonde_int{i}(1:length(rm));
            WV_sonde(:,i) = WV_sonde_int{i}(1:length(rm));
        end
%         rm_sonde(:,i) = rm_sonde_int{i};
%         T_sonde(:,i) = T_sonde_int{i};
%         P_sonde(:,i) = P_sonde_int{i};

        %interp sonde time
        [rm_sgp{i},IA,IC] = unique(rm_sgp{i});
        sonde_time(1:length(rm_sonde_int{i}),i) = interp1(rm_sgp{i},sondeStruc(i).time(IA),rm_sonde_int{i})';
        
        if length(sonde_time) < length(rm)
            sonde_time = [sonde_time; sonde_time(end).*ones(length(rm)-length(sonde_time),1)];
        end

         %Find index of sonde in time vector
         for j = 1:length((rm))
            [~, sonde_ind(j,i)]=min(abs(sonde_datetime(i)+seconds(sonde_time(j,i))-date_ts));
         end
    else
        sonde_ind = [];
    end
    
end



%%

% O2 Counts
o2on_intp = DataStructure2.MCS.Channel2.Data';
o2off_intp = DataStructure2.MCS.Channel10.Data';
o2on_intp_mol = DataStructure2.MCS.Channel0.Data';
o2off_intp_mol = DataStructure2.MCS.Channel8.Data';

% Count dead time correction
nsPerBin = 250e-9; %[s] nano seconds per bin
summedBins = 14000/2; %[] number of bins summed for one profile
% Deadtime for apd SPCM-780-12-FC
deadTime = 22e-9; %[s] 
o2on_countRate = o2on_intp /(nsPerBin * summedBins);
o2on_correctionFactor = 1./(1 - deadTime*o2on_countRate);
o2on_correctionFactor(o2on_correctionFactor<=0)=10;
o2on_intp = o2on_correctionFactor .* o2on_intp;
o2off_countRate = o2off_intp /(nsPerBin * summedBins);
o2off_correctionFactor = 1./(1 - deadTime*o2off_countRate);
o2off_correctionFactor(o2off_correctionFactor<=0)=10;
o2off_intp = o2off_correctionFactor .* o2off_intp;

o2on_countRate_mol = o2on_intp_mol /(nsPerBin * summedBins);
o2on_correctionFactor_mol = 1./(1 - deadTime*o2on_countRate_mol);
o2on_correctionFactor_mol(o2on_correctionFactor_mol<=0)=10;
o2on_intp_mol = o2on_correctionFactor_mol .* o2on_intp_mol;
o2off_countRate_mol = o2off_intp_mol /(nsPerBin * summedBins);
o2off_correctionFactor_mol = 1./(1 - deadTime*o2off_countRate_mol);
o2off_correctionFactor_mol(o2off_correctionFactor_mol<=0)=10;
o2off_intp_mol = o2off_correctionFactor_mol .* o2off_intp_mol;


% ======================
% === O2 Backscatter ===
% ======================
%Combined APD
A_co = 12.25;
RC_co = 0.1535;
%molecular APD
A_mol = 3.112;
RC_mol = 0.6873;


% --- O2 Background Subtraction ---
% Online
bg_o2on = mean(o2on_intp(end-20:end,:));% Take mean of last data points
o2on_bgsub = o2on_intp - bg_o2on;       % Background subtracted
o2on_bgsub(o2on_bgsub < 0) = 0;         % Minimum of zero
%o2on_noise = o2on_bgsub(1:i_range,:);   % Shorten to length of range vector
%o2on_bgsub = o2on_bgsub - (o2on_bgsub(3,:)-10).*exp(-rm_raw_o2/1000/RC_co); %Subtract recovery time model

bg_o2on_mol = mean(o2on_intp_mol(end-20:end,:));% Take mean of last data points
o2on_bgsub_mol = o2on_intp_mol - bg_o2on_mol;       % Background subtracted
o2on_bgsub_mol(o2on_bgsub_mol < 0) = 0;         % Minimum of zero
%o2on_noise_mol = o2on_bgsub_mol(1:i_range,:);   % Shorten to length of range vector
%o2on_bgsub_mol = o2on_bgsub_mol - (o2on_bgsub_mol(3,:)-10).*exp(-rm_raw_o2/1000/RC_mol);

% Offline
bg_o2off = mean(o2off_intp(end-20:end,:));% Take mean of last data points
o2off_bgsub = o2off_intp - bg_o2off;      % Background subtracted
o2off_bgsub(o2off_bgsub < 0) = 0;         % Minimum of zero
%o2off_noise = o2off_bgsub(1:i_range,:);   % Shorten to length of range vector
%o2off_bgsub = o2off_bgsub - (o2off_bgsub(3,:)-10).*exp(-rm_raw_o2/1000/RC_co); %Subtract recovery time model

bg_o2off_mol = mean(o2off_intp_mol(end-20:end,:));% Take mean of last data points
o2off_bgsub_mol = o2off_intp_mol - bg_o2off_mol;      % Background subtracted
o2off_bgsub_mol(o2off_bgsub_mol < 0) = 0;         % Minimum of zero
%o2off_noise_mol = o2off_bgsub_mol(1:i_range,:);   % Shorten to length of range vector



%%

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
t_avg = 31;                     %[min]

% Range oversample
oversample = 9;                 %[bins] Oversample must be even

oversample = 7;
t_avg = 45;

oversample = 7;
%oversample = 3;
t_avg = 21;

oversample = 1;
t_avg = 31;

% oversample = 3;
% t_avg = 1;

% Moving average in time and range
% Replicate edges to prevent problems caused by zero padding
% Even kernel causes data to shift by half step -- interpolate back to original time & range vectors
k = ones(oversample,t_avg)./(oversample*t_avg);     % Kernel
% k = triang(oversample).*k;
% k = k.*triang(t_avg)';
% k = k/trapz(trapz(k));
%kaiser window
%k = kaiser(oversample,0.5).*k;
%k = k.*kaiser(t_avg,0.5)';
%k = k/trapz(trapz(k));


% load('overlapFile2.mat')
% load('fullOverlapMolchannel.mat')
% 
% % range corrected returns
% o2on_noise = log(o2on_noise .* rm.^2  ./ overlapFile)./ Over;
% o2on_noise(o2on_noise<0) = 0;
% 
% o2off_noise = log(o2off_noise .* rm.^2  ./ overlapFile)./ Over;
% o2off_noise(o2off_noise<0) = 0;
% 
% o2on_noise_mol = log(o2on_noise_mol .* rm.^2 )./ Over;
% o2on_noise_mol(o2on_noise_mol<0) = 0;
% 
% o2off_noise_mol = log(o2off_noise_mol .* rm.^2 )./ Over;
% o2off_noise_mol(o2off_noise_mol<0) = 0;

% Online
% % o2on_noise_pad = padarray(o2on_noise,[oversample/2,t_avg/2],'replicate');
% % o2on_filt = filter2(k,o2on_noise_pad,'valid');
% % %o2on_filt = o2on_noise;
% % %o2on = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt,ts,rm);
% % o2on = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt(1:end-1,1:end-1),ts,rm);
o2on = filter2(k,o2on_noise,'same');
o2on = fillmissing(o2on,'nearest',1); % Fill in NaNs in dimension 1
o2on = fillmissing(o2on,'nearest',2); % Fill in NaNs in dimension 2

% % o2on_noise_pad_mol = padarray(o2on_noise_mol,[oversample/2,t_avg/2],'replicate');
% % o2on_filt_mol = filter2(k,o2on_noise_pad_mol,'valid');
% % %o2on_filt_mol = o2on_noise_mol;
% % %o2on_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt_mol,ts,rm);
% % o2on_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2on_filt_mol(1:end-1,1:end-1),ts,rm);
o2on_mol = filter2(k,o2on_noise_mol,'same');
o2on_mol = fillmissing(o2on_mol,'nearest',1); % Fill in NaNs in dimension 1
o2on_mol = fillmissing(o2on_mol,'nearest',2); % Fill in NaNs in dimension 2

% Offline
% % o2off_noise_pad = padarray(o2off_noise,[oversample/2,t_avg/2],'replicate');
% % o2off_filt = filter2(k,o2off_noise_pad,'valid');
% % %o2off_filt = o2off_noise;
% % %o2off = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt,ts,rm);
% % o2off = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt(1:end-1,1:end-1),ts,rm);
o2off = filter2(k,o2off_noise,'same');
o2off = fillmissing(o2off,'nearest',1); % Fill in NaNs in dimension 1
o2off = fillmissing(o2off,'nearest',2); % Fill in NaNs in dimension 2

% % o2off_noise_pad_mol = padarray(o2off_noise_mol,[oversample/2,t_avg/2],'replicate');
% % o2off_filt_mol = filter2(k,o2off_noise_pad_mol,'valid');
% % %o2off_filt_mol = o2off_noise_mol;
% % %o2off_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt_mol,ts,rm);
% % o2off_mol = interp2(ts-t_step/2,rm-rangeBin/2,o2off_filt_mol(1:end-1,1:end-1),ts,rm);
o2off_mol = filter2(k,o2off_noise_mol,'same');
o2off_mol = fillmissing(o2off_mol,'nearest',1); % Fill in NaNs in dimension 1
o2off_mol = fillmissing(o2off_mol,'nearest',2); % Fill in NaNs in dimension 2

%o2on = o2on_mol;
%o2off = o2off_mol;



% %unfiltered
% o2on = o2on_noise;
% o2off = o2off_noise;
% o2on_mol = o2on_noise_mol;
% o2off_mol = o2off_noise_mol;

%convert to log rang
% o2off_noise = log(o2off_noise.*rm.^2);
% o2on_noise = log(o2on_noise.*rm.^2);
% o2off_noise_mol = log(o2off_noise_mol.*rm.^2);
% o2on_noise_mol = log(o2on_noise_mol.*rm.^2);

%sgolay filtering test
% % % order = 0;
% % % o2off_sg = sgolayfilt(o2off_noise, order, t_avg+1,[],2);
% % % o2off = sgolayfilt(o2off_sg, order, oversample+1,[],1);
% % % o2off(o2off<0)=0;
% % % 
% % % o2off_sg_mol = sgolayfilt(o2off_noise_mol, order, t_avg+1,[],2);
% % % o2off_mol = sgolayfilt(o2off_sg_mol, order, oversample+1,[],1);
% % % o2off_mol(o2off_mol<0)=0;
% % % o2on_sg = sgolayfilt(o2on_noise, order, t_avg+1,[],2);
% % % o2on = sgolayfilt(o2on_sg, order, oversample+1,[],1);
% % % o2on(o2on<0)=0;
% % % o2on_sg_mol = sgolayfilt(o2on_noise_mol, order, t_avg+1,[],2);
% % % o2on_mol = sgolayfilt(o2on_sg_mol, order, oversample+1,[],1);
% % % o2on_mol(o2on_mol<0)=0;


% o2off = exp(o2off)./rm.^2;
% o2on = exp(o2on)./rm.^2;
% o2off_mol = exp(o2off_mol)./rm.^2;
% o2on_mol = exp(o2on_mol)./rm.^2;
% 
% o2off_noise = exp(o2off_noise)./rm.^2;
% o2on_noise = exp(o2on_noise)./rm.^2;
% o2off_noise_mol = exp(o2off_noise_mol)./rm.^2;
% o2on_noise_mol = exp(o2on_noise_mol)./rm.^2;


%%
% Apply online overlap correction
% load('overlapFile2.mat')
% o2on_mol = o2on_mol .* overlapFile;
% o2on = o2on + o2on_mol;

% %Apply online overlap correction
% % % load('overlapFile3.mat')
% % % o2on_mol = o2on_mol .* overlapFile;
%o2on = o2on + o2on_mol;

%Apply online overlap correction
% load('overlapFile4.mat')
% o2on_mol = o2on_mol .* overlapFile;
% o2on = o2on + o2on_mol;


%%

%======================
%= Cloud and SNR mask =
%======================
disp('Calculating masks')
cloud_p_point = 18;
SNR_threshold = 2;
SD_threshold = 6*10^9;
SD_threshold = 10*10^8;
SD_threshold = 1*10^9;
SD_threshold = 5*10^8;
%SD_threshold = 1*10^8;
%SD_threshold = 5*10^12;
[SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2(o2on,o2off,rm,ts,cloud_p_point,SNR_threshold,SD_threshold,oversample,t_avg);

%%
disp('Calculating model')

% Calculating temperature and pressure model

Ts = weather_Temperature_interp + 273.15 ;          %surface temperature from weather station [K]
Ps = weather_absPressure_interp / 1013.25;         %absolute surface pressure from weather station [atm]


lapseRate = -6.5;                                   %[K/km] Guess adiabatic lapse rate  typically -6.5 up to 10km
%lapseRate = -15;
%lapseRate = -9;
lapseRate = lapseRate / 1000;                       %[K/m] 

T = Ts + lapseRate .* rm;                           %[K] (1 x r) Temperature model as a function of r 

P = Ps .* (Ts./T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  

%water vapor
Pws = exp(77.3450+0.0057.*T-7235./T)./T.^8.2;
%WV = weather_WV_interp.*0.0022.*Pws./T./100;
%WV = weather_WV_interp.*0.0022.*Pws./T./100;
WV = weather_WV_interp.*Pws./T./100/kb;%water vapror molc/m^3


%%%%%%%WV = ones(i_range,i_time).*WV_sonde(:,1);


%WV=WV;%.*50;

%WV(:,:)=0;
%%

% for ii = 1:2
%     fprintf('Large Loop %d\n',ii)
%     if ii>1
%         for iii =1:i_range
%             for jjj =1:i_time
%                 if ~isnan(T_finalm(iii,jjj))
%                     T(iii,jjj) = T_finalm(iii,jjj);
%                 end
%             end
%         end
%     end
    
   
%%
% Use ncar reanalysis data
%bozeman lat and longitude
% % % lat = 45.7;
% % % lon = -111;
% % % bzElevation = 1461;%[m]
% % % [Tnoaa,rmnoaa,latRe,lonRe,hgtReanalysys] = tempReanalysis(Ts,Ps,path,span_days(1),ts,lat,lon);
% % % rmnoaa = rmnoaa - bzElevation;
% for i =1:length(ts)
%     rmnoaashort(:,i) = rmnoaa(rmnoaa(:,i)<=7100 & rmnoaa(:,i)>=70,i);
%     Tnoaashort(:,i) = Tnoaa(rmnoaa(:,i)<=7100 & rmnoaa(:,i)>=70,i);
% end

%T = interp1(rmnoaashort(:,1),Tnoaashort,rm);
%T = fillmissing(T,'nearest');
%P = Ps .* (Ts./T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
%%
%-----Create Spectrum vectors----

%lambda_online = 769.2330;                          %[nm] Cobleigh online wavelength
%lambda_offline = 769.3184;                         %[nm] Cobleigh offline wavelength
% lambda_online = O2Online_Wavelength'*10^9;          %[nm] Updated online wavelength
% lambda_offline = O2Offline_Wavelength'*10^9;        %[nm] Updated offline wavelength
lambda_online = interp1(Options.TimeGrid,DataStructure2.Laser.O2Online.WavelengthActual,thr);
lambda_offline = interp1(Options.TimeGrid,DataStructure2.Laser.O2Offline.WavelengthActual,thr);
%shorten to length of time vector if needed
%lambda_online = lambda_online(1:length(DataStructure2.MCS.Channel2.Data'));
%lambda_offline = lambda_offline(1:length(DataStructure2.MCS.Channel2.Data'));
lambda_online = lambda_online(1:length(ts));
lambda_offline = lambda_offline(1:length(ts));

lambda_online(1:length(ts)) = 769.7958;
lambda_offline(1:length(ts)) = 770.1085;

%lambda_online(1:length(ts)) = 769.7954;
%lambda_offline(1:length(ts)) = 770.1081;
nu_online = 10^7./lambda_online;                    %[cm-1] Online wavenumber
nu_offline = 10^7./lambda_offline;                  %[cm-1] Offline wavenumber

nu01 = nu_online;                                   %[cm-1] Set center of scan to 
nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
nuBin = 0.00222;                                    %[cm-1] Scan increment
%nuBin = nuBin/2;
%nuBin = 8.8800e-04;
%nuBin = 0.0007;                                    %[cm-1] Scan increment
nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nu01_off = nu_offline;
nuMin_off = nu01_off-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = nu01_off+0.334;                                 %[cm-1] Scan upper bound
nu_scan_off = (nuMin_off:nuBin:nuMax_off);

% %new nu
% lambdaCenter = 769.7958;
% lambdaWidth = 0.2;
% lambdaDelta = .00002;
% lambda = lambdaCenter-lambdaWidth/2:lambdaDelta:lambdaCenter+lambdaWidth/2;%[nm]
% nu_scan = 10^7./lambda;
% %nu_scan = permute(nu,[3 2 1]);


%i_scan = length(nu_scan);                           %[none] length of scan vector

lambda_scan = 10^7./nu_scan;                        %[nm](1 x lambda) Scan vector
%f_scan = nu_scan * c * 100;                         %[Hz]( x f) Scan vector

nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension

lambda_scan_3D_short = 10^7./nu_scan_3D_short;
lambda_scan_3D_short_off = 10^7./nu_scan_3D_short_off;
i_scan_3D_short = length(nu_scan_3D_short);         %[none] length of scan vector

del_nu = nu_scan_3D_short-nu_online;                %[1/cm] difference from center
del_lambda = lambda_scan_3D_short-lambda_online;
%del_f = del_nu*100*c*10^-9;                         %[GHz] difference from center


[~,online_index] = min(abs(nu_online - nu_scan_3D_short),[],3);%finding index of online wavenumber

[~,offline_index] = min(abs(nu_offline - nu_scan_3D_short_off),[],3);%finding index of online wavenumber



 
absorption = absorption_O2_770_model(T,P,nu_online,WV); %[m-1] Funcrtion to calculate theoretical absorption
%absorption_off = absorption_O2_770_model(T,P,nu_offline);

% for i=1:numel(sonde_datetime)
%     for j = 1:length(rm)
%         if isdatetime(sonde_datetime(i))
%             absorption_sonde{i}(j) = absorption_O2_770_model(T_sonde(j,i),P_sonde(j,i),nu_online(sonde_ind(j,i)),WV(j,sonde_ind(j,i))); %[m-1] Funcrtion to calculate theoretical absorption
%         else
%             absorption_sonde{i}(j) = nan(i_range,1);
%             T_sonde = nan(i_range,1);
%             P_sonde = nan(i_range,1);
%         end
%     end
% end

for i=1:numel(sonde_datetime)
        if isdatetime(sonde_datetime(i))
            absorption_sonde{i} = diag(absorption_O2_770_model(T_sonde(:,i),P_sonde(:,i),nu_online(sonde_ind(:,i)),WV(:,sonde_ind(:,i)))); %[m-1] Funcrtion to calculate theoretical absorption
        else
            absorption_sonde{i} = nan(i_range,1);
            T_sonde = nan(i_range,1);
            P_sonde = nan(i_range,1);
        end
end

%%
%=====================
%= Backscatter ratio =
%=====================

%[Bm,Ba,BR,Te,Pe] = BackscatterRatio(ts',0,rm'/1000,0,o2off',o2off_mol',weather_absPressure_interp',weather_Temperature_interp',lambda_offline');
%load('overlapcorrection0619.mat')


%[Bm,Ba,BR]= BackscatterRatiov2(ts,rm,o2off,o2off_mol,T,P,lambda_offline);
%[Bm,Ba,BR]= BackscatterRatioV3(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
%[Bm,Ba,BR]= BackscatterRatioV3(ts,rm,o2off,o2off_mol,T,P,lambda_offline);


if span_days(1)<datetime(2020,10,6,'TimeZone','UTC')
    load('overlap.mat')
    overlapcorrection = interp1(rm_raw_o2,Correction,rm);
    [Bm,Ba,BR]= BackscatterRatioV3(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
    
    %BR = interp2(ts,rm+50,BR,ts,rm,'nearest',1);
    
        elseif span_days(1)>=datetime(2021,3,12,'TimeZone','UTC')
    %addpath '\BSR retrieval'
    %load('Overlap0304.mat')
    %LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
    LidarData.Range = rm;
    LidarData.Time = ts;
    LidarData.OfflineCombinedAverageCounts = o2off;
    LidarData.OfflineMolecularAverageCounts = o2off_mol;
    WeatherData.Temperature = T;
    WeatherData.Pressure = P;
    [LidarData]=BackscatterRetrievalRayleighBrillouin0312(LidarData,WeatherData);
    BR = LidarData.BackscatterRatio;
    
    
    elseif span_days(1)>datetime(2021,2,20,'TimeZone','UTC')
    %addpath '\BSR retrieval'
    %load('Overlap0304.mat')
    %LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
    LidarData.Range = rm;
    LidarData.Time = ts;
    LidarData.OfflineCombinedAverageCounts = o2off;
    LidarData.OfflineMolecularAverageCounts = o2off_mol;
    WeatherData.Temperature = T;
    WeatherData.Pressure = P;
    [LidarData]=BackscatterRetrievalRayleighBrillouin(LidarData,WeatherData);
    BR = LidarData.BackscatterRatio;
    
elseif span_days(1)>datetime(2021,1,28,'TimeZone','UTC')
        load('Overlap1104.mat')
    overlapcorrection = interp1(rm_raw_o2,overlapcorrection,rm);
    %[Bm,Ba,BR]= BackscatterRetrievalRayleigh(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
    LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
    LidarData.OfflineMolecularAverageCounts = o2off_mol;
    WeatherData.Temperature = T;
    WeatherData.Pressure = P;
    [LidarData]=BackscatterRetrieval_2_10_21(LidarData,WeatherData);
    BR = LidarData.BackscatterRatio;
        


else
    load('Overlap1006.mat')
    overlapcorrection = interp1(rm_raw_o2,overlapcorrection,rm);
    [Bm,Ba,BR]= BackscatterRatioV4(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);

    load('Overlap1104.mat')
    overlapcorrection = interp1(rm_raw_o2,overlapcorrection,rm);
    %[Bm,Ba,BR]= BackscatterRetrievalRayleigh(ts,rm,o2off./overlapcorrection,o2off_mol,T,P,lambda_offline);
    LidarData.OfflineCombinedAverageCounts = o2off./overlapcorrection;
    LidarData.OfflineMolecularAverageCounts = o2off_mol;
    WeatherData.Temperature = T;
    WeatherData.Pressure = P;
    [LidarData]=BackscatterRetrieval(LidarData,WeatherData);
    BR = LidarData.BackscatterRatio;
end

% Averaging
% % BR_noise = padarray(BR,[oversample/2,t_avg/2],'replicate');
% % BR_filt = filter2(k,BR_noise,'valid');
% % BR = interp2(ts-t_step/2,rm-rangeBin/2,BR_filt(1:end-1,1:end-1),ts,rm);
% % BR = fillmissing(BR,'nearest',1); % Fill in NaNs in dimension 1
% % BR = fillmissing(BR,'nearest',2); % Fill in NaNs in dimension 2

%BR=smoothdata(BR,1,'g',30);


%BSR = interp2(ts,rm,BR,ts,rm-3*rangeBin,'spline');
BSR = BR;

%BSR = 1.15*BSR;

%make upper BSR=1
% % BSRmeanEnd=mean(BSR(end-20,end),1);
% % BSR = BSR./BSRmeanEnd;

% cloud_SDm_above = masking(o2off',o2off_mol',BR);
% cloud_SDm_above = cloud_SDm_above';

%%

%==========================
%= Pertabative absorption =
%==========================

disp('Calculating absorption')

%Masking couts
% o2onm = o2on .* cloud_SDm_above;
% o2onm(o2onm<=0) = NaN;
% for i = 1:i_time
%     if ismissing(o2onm(end,i))
%         o2onm(end,i)=0;
%     end
% end
% o2onm = fillmissing(o2onm,'linear');
% 
% o2offm = o2off .* cloud_SDm_above;
% o2offm(o2offm<=0) = NaN;
% for i = 1:i_time
%     if ismissing(o2offm(end,i))
%         o2offm(end,i)=0;
%     end
% end
% o2offm = fillmissing(o2offm,'linear');

% === Zeroth Order ===
 ind_r_lo = 1:i_range-oversample;                                            % High range vector
 ind_r_hi = 1+oversample:i_range;                                            % Low range vector
% % % ind_r_lo = 1:i_range-1;                                            % High range vector
% % % ind_r_hi = 1+1:i_range;                                            % Low range vector
 ln_o2 = log((o2on(ind_r_lo,:).*o2off(ind_r_hi,:))./(o2on(ind_r_hi,:).*o2off(ind_r_lo,:))); % Natural log of counts
%ln_o2 = log((o2onm(ind_r_lo,:).*o2offm(ind_r_hi,:))./(o2onm(ind_r_hi,:).*o2offm(ind_r_lo,:))); % Natural log of counts
%ln_o2_total = interp1(rm(ind_r_lo),ln_o2,rm,'nearest','extrap'); % For plotting ln_o2 vs rm

% Zeroth order term
 alpha_0_raw = ln_o2./2./(rangeBin*oversample);                              %[1/m] 
% % % alpha_0_raw = ln_o2./2./(rangeBin);                              %[1/m] 

alpha_0 = interp2(ts,rm(ind_r_lo),alpha_0_raw,ts,rm);

% --- Filter alpha_0 ---
% % Savitzky-Golay filtering
% alpha_0_sgfilt = sgolayfilt(alpha_0_raw, 1, 2*oversample+1);

% Moving average
% % % alpha_0_pad = padarray(alpha_0_raw,[oversample/2,t_avg/2],'replicate');
% % % alpha_0_filt = filter2(k,alpha_0_pad,'valid');
% % % %%%alpha_0 = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2-rangeBin*oversample/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% % % alpha_0 = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% % % %alpha_0 = interp2(ts,rm(ind_r_lo),alpha_0_raw,ts,rm);
% % % alpha_0 = fillmissing(alpha_0,'nearest',1);                                 % Fill in NaNs in dimension 1
% % % alpha_0 = fillmissing(alpha_0,'nearest',2);                                 % Fill in NaNs in dimension 2
% % % %alpha_0(alpha_0==Inf | alpha_0==-Inf) = 0;


% % % %     disp('SG derivative')
% % % %     int_der = -log(o2on)+log(o2off);
% % % %     tic
% % % %     [b,g] = sgolay(2,oversample);
% % % %     parfor j=1:i_time
% % % %         alpha_0(:,j) = conv(int_der(:,j), factorial(1)/(-rangeBin)^1 * g(:,2), 'same')/2;
% % % %     end
% % % %     toc
% % % %     
% % % %     %alpha 0 mol
% % % % 
% % % %     int_der = -log(o2on_mol)+log(o2off_mol);
% % % %     tic
% % % %     [b,g] = sgolay(2,oversample);
% % % %     parfor j=1:i_time
% % % %         alpha_0_mol(:,j) = conv(int_der(:,j), factorial(1)/(-rangeBin)^1 * g(:,2), 'same')/2;
% % % %     end
% % % %     
% % % %  
% % % %     o2off_mol_corr =o2off_mol.*BSR;
% % % %     int_der = -log(o2on_mol)+log(o2off_mol_corr);
% % % %     tic
% % % %     [b,g] = sgolay(2,oversample);
% % % %     parfor j=1:i_time
% % % %         alpha_0_mol_corr(:,j) = conv(int_der(:,j), factorial(1)/(-rangeBin)^1 * g(:,2), 'same')/2;
% % % %     end
% % % %     alpha_0_mol_corr = real(alpha_0_mol_corr);
    
    %alpha_0 = (alpha_0+alpha_0_mol_corr)/2;
    %alpha_0 = (.8*alpha_0+0.2*alpha_0_mol_corr);


%%
% ind_r_lo = 1:i_range-1;                                            % High range vector
% ind_r_hi = 1+1:i_range;                                            % Low range vector
% ln_o2 = log((o2on(ind_r_lo,:).*o2off(ind_r_hi,:))./(o2on(ind_r_hi,:).*o2off(ind_r_lo,:))); % Natural log of counts
% %ln_o2_total = interp1(rm(ind_r_lo),ln_o2,rm,'nearest','extrap'); % For plotting ln_o2 vs rm
% 
% % Zeroth order term
% alpha_0_raw = ln_o2./2./(rangeBin*1);                              %[1/m] 
% 
% % --- Filter alpha_0 ---
% % % Savitzky-Golay filtering
% % alpha_0_sgfilt = sgolayfilt(alpha_0_raw, 1, 2*oversample+1);
% 
% % Moving average
% alpha_0_pad = padarray(alpha_0_raw,[oversample/2,t_avg/2],'replicate');
% alpha_0_filt = filter2(k,alpha_0_pad,'valid');
% alpha_0 = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% %alpha_0 = interp2(ts,rm(ind_r_lo),alpha_0_raw,ts,rm);
% alpha_0 = fillmissing(alpha_0,'nearest',1);                                 % Fill in NaNs in dimension 1
% alpha_0 = fillmissing(alpha_0,'nearest',2);                                 % Fill in NaNs in dimension 2


% % ------ fitting alpha zero ---------
% alpha_0_fit = zeros(i_range,i_time);
%  for i = 1:i_time
%      alpha_0_fit_obj = fit(rm(26:106,1),alpha_0(26:106,i),'poly1');
%      alpha_0_fit(:,i) = alpha_0_fit_obj(rm);
%  end

%--first order--
 
% fsr_O2 = 157.9;%[GHz] etalon free spectral range
% finesse_O2 = 15.43;%etalon finesse
% %filter_bw = 13;%[nm]bandwidth of filter
% %filter_bw2 = 1;%[nm]bandwidth of filter
% %NBF_blocking = 10^-6;%NBF out of band blocking
% %NBF_blocking2 = 10^-4;%NBS out of band blocking
% %NBF_transmission = 0.9;
% %NBF_transmission2 = 0.7;
% 
% % ======================================
% % === Normalized Etalon Transmission ===
% % ======================================
% % --- O2 channel ---
% fsr_O2_nm = lambda_online.^2 * fsr_O2 / c;                          %[nm] Free spectral range converted to nm
% delta_phase = 2 * pi * lambda_online.^2 ./ fsr_O2_nm ./ lambda_scan_3D_short; %[radians] Phase difference between orders
% delta_phase = delta_phase - 2 * pi * lambda_online / fsr_O2_nm;     %[radians] Shifting phase so that zero is at the online wavelength
% F_O2 = 1 ./ (sin(pi/2/finesse_O2)^2);                               %[none] Coefficient of finesse
% T_etalon = 1 ./ (1 + F_O2 * sin(delta_phase/2).^2);                 %[none] Etalon transmission as a function of wavelength

%load('etalonTransmission.mat')
%load('wavelengthEtalonScan.mat')
% load('etalonTransmission2.mat')
% load('wavelengthEtalonScan2.mat')
% %[wavelengthEtalonScan,IA,IC] = unique(wavelengthSorted);
% [wavelengthEtalonScan,IA,IC] = unique(wavelengthEtalonScan);
% etalonTransmission = etalonTransmission(IA);
% T_etalon = interp1(wavelengthEtalonScan,etalonTransmission,lambda_scan_3D_short);

%%
load('TransmittanceData.mat')
%load('EtalonScan.mat')
T_etalon_on = double(interp1(double(OnlineWavelength)*10^9,OnlineCombinedTransmittance,lambda_scan_3D_short));
T_etalon_off = double(interp1(double(OfflineWavelength)*10^9,OfflineCombinedTransmittance,lambda_scan_3D_short_off));
%%

% % % --- Spectral distribution using the initial temperature profile guess ---
% % % c_doppler_O2 = m_air*c^2./(8*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
% % % doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T); %[m] Doppler broadended lineshape      
% % % c_doppler_O2_off = m_air*c^2./(8*(nu_offline*100).^2*kb);                   %[m^2 K] Doppler coefficient
% % % doppler_O2_un_ret_off = ((c_doppler_O2_off./T/pi).^0.5).*exp(-c_doppler_O2_off.*(nu_offline*100-nu_scan_3D_short_off*100).^2./T); %[m] Doppler broadended lineshape    
% % 
% % cB = 1.2;%Brullouion correction to doppler gaussian half width
% % %cB = 1;%Brullouion correction to doppler gaussian half width
% % 
% % cB = -0.01*(rkm+1.5) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km
% % %cB=1;
% % 
% % c_doppler_O2 = m_air*c^2./(8*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
% % doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         
% % c_doppler_O2_off = m_air*c^2./(8*(nu_offline*100).^2*kb);                   %[m^2 K] Doppler coefficient
% % doppler_O2_un_ret_off = ((c_doppler_O2_off./T/pi).^0.5).*exp(-c_doppler_O2_off.*(nu_offline*100-nu_scan_3D_short_off*100).^2./T./cB.^2); %[m] Doppler broadended lineshape    
% % 
% % 
% % 
% % % c_doppler_O2 = m_air*c^2./(8*pi*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
% % % doppler_O2_un_ret = ((c_doppler_O2./T).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T); %[m] Doppler broadended lineshape     
% % norm_O2_ret = trapz(doppler_O2_un_ret,3).*nuBin*100;                   %[none] Lineshape integral
% % doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape
% % 
% % norm_O2_ret_off = trapz(doppler_O2_un_ret_off,3).*nuBin*100;                   %[none] Lineshape integral
% % doppler_O2_ret_off = doppler_O2_un_ret_off./norm_O2_ret;                       %[m] Normalized doppler lineshape
% % 
% % % Check if doppler_o2_ret is normalized to 1 when integrated across frequency
% % doppler_o2_ret_check = trapz(doppler_O2_ret,3).*nuBin*100;              %[none]
% % doppler_o2_ret_check_off = trapz(doppler_O2_ret_off,3).*nuBin*100;
% % 
% % % --- Backscatter Lineshape g ---
% % g1_m_on = 1./BSR .* doppler_O2_ret;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
% % g1_m_off = 1./BSR .* doppler_O2_ret_off;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
% % g1_a_on = zeros(i_range,i_time,i_scan_3D_short);                       % Initalize aerosol lineshape
% % g1_a_off = zeros(i_range,i_time,i_scan_3D_short);                       % Initalize aerosol lineshape
% % for i = 1:i_time
% %     g1_a_on(:,i,online_index(i)) = (1 - 1./BSR(:,i)) / nuBin / 100 ; %[m] aerosol backscatter lineshape
% %     g1_a_off(:,i,offline_index(i)) = (1 - 1./BSR(:,i)) / nuBin / 100 ; %[m] aerosol backscatter lineshape
% % end
% % g1 = g1_a_on + g1_m_on;                                                   %[m] Combined backscatter lineshape
% % g1_off = g1_a_off + g1_m_off;                                                   %[m] Combined backscatter lineshape
% % 
% % g1_check = trapz(g1,3).*nuBin*100;                                %[none] Check if integral of g1 is normalized to 1
% % g1_check_off = trapz(g1_off,3).*nuBin*100;                                %[none] Check if integral of g1 is normalized to 1
% % 
% % %derivative of lineshape dg/dr
% % % ind_r_lo = 1:i_range-1;                                            % High range vector
% % % ind_r_hi = 1+1:i_range;                                            % Low range vector
% % dg1_dr = 1*(g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(rangeBin*oversample); %[none] Derivative over oversamped range
% % % % dg1_dr = (g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(rangeBin); %[none] Derivative over oversamped range
% % dg1_dr = interp1(rm(ind_r_lo),dg1_dr,rm,'nearest','extrap');         %[none] Make dg/dr the same size as g
% % 
% % dg1_dr_off = 1*(g1_off(ind_r_hi,:,:) - g1_off(ind_r_lo,:,:)) ./(rangeBin*oversample); %[none] Derivative over oversamped range
% % % % dg1_dr_off = (g1_off(ind_r_hi,:,:) - g1_off(ind_r_lo,:,:)) ./(rangeBin); %[none] Derivative over oversamped range
% % dg1_dr_off = interp1(rm(ind_r_lo),dg1_dr_off,rm,'nearest','extrap');         %[none] Make dg/dr the same size as g
% % %%
% % tic
% % disp('calculation absorption f')
% % absorption_f = absorption_O2_770_model_wavenumber(T,P,nu_scan_3D_short,WV); %[m] lineshape function 
% % toc
% % %%
% % % tic
% % % disp('calculation absorption f2')
% % % [absorption_f2,cross_section] = absorption_O2_770_PCA(T,P,nu_scan_3D_short,WV);
% % % toc
% % %%
% % tic
% % disp('calculating f')
% % f = ones(size(absorption_f));
% % for i = 1:i_time
% %     %f(:,i,:) = absorption_f(:,i,:) ./ absorption_f(:,i,online_index(i));                  %[none] Normalize cross section to line center
% %     f(:,i,:) = absorption_f(:,i,:) ./ max(absorption_f(:,i,:),[],3);    %[none] Normalize cross section to line center
% % end
% % toc
% % %%
% % % tic
% % % disp('calculating f2')
% % % for i = 1:i_time
% % %     %f(:,i,:) = absorption_f(:,i,:) ./ absorption_f(:,i,online_index(i));                  %[none] Normalize cross section to line center
% % %     f2(:,i,:) = absorption_f2(:,i,:) ./ max(absorption_f2(:,i,:),[],3);    %[none] Normalize cross section to line center
% % % end
% % % toc
 %%    
% % % --- Zeroth Order Transmission ---
% % Tm0 = exp(-cumtrapz(rm,alpha_0.*f,1));      %[none] Zeroth order transmission
% % 
% % % Integrand terms
% % % Online
% % zeta = g1.*T_etalon_on;                        %[m]
% % eta = dg1_dr.*T_etalon_on;                     %[none]
% % 
% % % Integrated terms
% % % Online
% % zeta_int = trapz(zeta.*Tm0,3)*nuBin*100;              %[none]
% % eta_int = trapz(eta.*Tm0,3)*nuBin*100;                %[1/m]
% % zeta_ls_int = trapz(zeta.*Tm0.*(1-f),3)*nuBin*100;    %[none]
% % % Offline
% % zeta_off = g1_off.*T_etalon_off;                        %[m]
% % eta_off = dg1_dr_off.*T_etalon_off;                     %[none]
% % zeta2_int = trapz(zeta_off,3)*nuBin*100;                  %[none]
% % eta2_int = trapz(eta_off,3)*nuBin*100;                    %[1/m]
% % 
% % 
% % % === First Order ===
% % W1 = zeta_ls_int./zeta_int;                 %[none]
% % G1 = eta_int./zeta_int - eta2_int./zeta2_int;%[1/m]
% % 
% % alpha_1_raw = 0.5.*(alpha_0_fit.*W1 + G1);  %[1/m]
% % alpha_1_raw = 0.5.*(alpha_0.*W1 + G1);      %[1/m]
% % 
% % % Moving average
% % % % % alpha_1_pad = padarray(alpha_1_raw,[oversample/2,t_avg/2],'replicate');
% % % % % alpha_1_filt = filter2(k,alpha_1_pad,'valid');
% % % % % alpha_1 = interp2(ts-t_step/2,rm-rangeBin/2,alpha_1_filt(1:end-1,1:end-1),ts,rm);
% % % % % alpha_1 = fillmissing(alpha_1,'nearest',1); % Fill in NaNs in dimension 1
% % % % % alpha_1 = fillmissing(alpha_1,'nearest',2); % Fill in NaNs in dimension 2
% %  alpha_1 = alpha_1_raw;
% % 
% % % alpha_1s = zeros(i_range,i_time);
% % % % ----- smoothing alpha 1 -----
% % % for i = 1:i_time
% % %     alpha_1s(:,i) = 1.15.*smooth(alpha_1(:,i), 4*oversample+1);
% % % end
% % 
% % % --- First Order Transmission Tm1 ---
% % Tm1 = exp(-cumtrapz(rm,oversample.*alpha_1.*f,1));      %[none] First order transmission
% % 
% % % === Second Order ===
% % % Integrated terms
% % zeta_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1),3)*nuBin*100;             %[none]
% % eta_Tm1_int = trapz(eta.*Tm0.*(1-Tm1),3)*nuBin*100;               %[1/m]
% % zeta_ls_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1).*(1-f),3)*nuBin*100;   %[none]
% % 
% % clear Tm0 Tm1
% % 
% % W2 = (zeta_ls_int.*zeta_ls_Tm1_int./(zeta_int.^2)) - (zeta_ls_Tm1_int./zeta_int);   %[none]
% % G2 = (eta_int.*zeta_Tm1_int./(zeta_int.^2)) - (eta_Tm1_int./zeta_int);              %[1/m]
% % 
% % alpha_2_raw = 0.5.*(alpha_1.*W1 + alpha_0_fit.*W2 + G2);
% % alpha_2_raw = 0.5.*(alpha_1.*W1 + alpha_0.*W2 + G2);    %[1/m]
% % 
% % %Moving average
% % % % % alpha_2_pad = padarray(alpha_2_raw,[oversample/2,t_avg/2],'replicate');
% % % % % alpha_2_filt = filter2(k,alpha_2_pad,'valid');
% % % % % alpha_2 = interp2(ts-t_step/2,rm-rangeBin/2,alpha_2_filt(1:end-1,1:end-1),ts,rm);
% % % % % alpha_2 = fillmissing(alpha_2,'nearest',1); % Fill in NaNs in dimension 1
% % % % % alpha_2 = fillmissing(alpha_2,'nearest',2); % Fill in NaNs in dimension 2
% %  alpha_2 = alpha_2_raw;
 
%  figure(666452)
%  pp_point = 840;
%  plot(alpha_0(:,pp_point),rm,alpha_0(:,pp_point)+alpha_1(:,pp_point),rm,alpha_0(:,pp_point)+alpha_1(:,pp_point)+alpha_2(:,pp_point),rm)
%  
 altitude = 1.5;%altitude in km
 %[alpha_1, alpha_2] = pertAbsorption(alpha_0, T_etalon_on, T, P, rm,rkm,m_air, nu_online, nu_scan_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,c,kb,altitude);
 [alpha_1, alpha_2] = pertAbsorption(alpha_0, T_etalon_on, T, P, rm,ts,rkm,m_air, nu_online, nu_scan_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,t_avg,c,kb,altitude);
%  figure(666453)
%  plot(alpha_0(:,pp_point),rm,alpha_0(:,pp_point)+alpha_1(:,pp_point),rm,alpha_0(:,pp_point)+alpha_1(:,pp_point)+alpha_2(:,pp_point),rm)
%  
 
% ----- smoothing alpha 2 -----
% alpha_2s = zeros(i_range,i_time);
% for i = 1:i_time
%     alpha_2s(:,i) = 1.15.*smooth(alpha_2(:,i), 4*oversample+1);
% end

% === Total alpha ===
alpha_total_raw = alpha_0 + alpha_1 + alpha_2;

% Force total alpha to its modeled surface value
%[~,cut] = min(abs(rm-1000));             % Index where rm is closest to chosen value
[~,cut] = min(abs(rm-500));             % Index where rm is closest to chosen value
[~,cut] = min(abs(rm-0));             % Index where rm is closest to chosen value
alpha_total = [absorption(1,:); NaN((cut - 1),i_time); alpha_total_raw(cut:end,:)];
alpha_total = fillmissing(alpha_total,'linear');
alpha_total = alpha_total(2:end,:);     % Remove surface value since rm starts at del_r

% Moving average
% % % alpha_total_pad = padarray(alpha_total,[oversample/2,t_avg/2],'replicate');
% % % alpha_total_filt = filter2(k,alpha_total_pad,'valid');
% % % alpha_total = interp2(ts-t_step/2,rm-rangeBin/2,alpha_total_filt(1:end-1,1:end-1),ts,rm);
% % % alpha_total = fillmissing(alpha_total,'nearest',1); % Fill in NaNs in dimension 1
% % % alpha_total = fillmissing(alpha_total,'nearest',2); % Fill in NaNs in dimension 2

%only do zeroth order
%alpha_total = alpha_0;

% ----- smoothing alpha total -----
% for i = 1:i_time
%     alpha_total_pad2(:,i) = padarray(alpha_total(:,i),[4*oversample+1,0],'replicate');
%     alpha_total_pad2(:,i) = smooth(alpha_total_pad2(:,i), 4*oversample+1);
%     alpha_total(:,i) = interp1(rm-rangeBin/2,alpha_total_pad2(4*oversample+1:end-(4*oversample+1)-1,i),rm);
% end

% for i = 1:i_time
%     alpha_total(:,i) = smooth(alpha_total(:,i), 4*oversample+1);
% end

%%
%apply SNR mask
alpha_0m = alpha_0 .* SNRm;
alpha_0m(alpha_0m == 0) = NaN;                  % Replace mask with NaNs

alpha_1 = alpha_1 .* SNRm;
alpha_1(alpha_1 == 0) = NaN;                    % Replace mask with NaNs

alpha_totalm = alpha_total .* cloud_SDm_above; %.* cloud_SDm_above;
%alpha_totalm(alpha_totalm <= 0) = NaN;          % Replace mask with NaNs

%alpha_total_raw = alpha_total_raw .* SNRm;
%alpha_total_raw(alpha_total_raw == 0) = NaN;    % Replace mask with NaNs

% % % % % for i = 1:i_time
% % % % %     alpha_totalm(:,i) = smooth(alpha_totalm(:,i), 4*oversample+1);
% % % % % end


% % for i = 1:i_time
% %     alpha_totalm(:,i) = smooth(alpha_totalm(:,i), 4*oversample+1);
% % end

%alpha_totalm = alpha_totalm .* SNRm .* cloud_SDm_above;
%alpha_totalm(alpha_totalm <= 0) = NaN;          % Replace mask with NaNs


%%
disp('Temp retrieval function')
%WV = 0; % set water vapor concentration to zero
tic
[T_final_test,L_fit_sm_test,Ts_fit,Patm_final,mean_lapse_rate,exclusion] =  temperatureRetrieval(T,ts,rm,P,WV,nu_online,alpha_totalm,SNRm,cloud_SDm_above);
toc



% % % alpha_total_pad = padarray(alpha_total,[oversample/2,t_avg/2],'replicate');
% % % alpha_total_filt = filter2(k,alpha_total_pad,'valid');
% % % alpha_total = interp2(ts-t_step/2,rm-rangeBin/2,alpha_total_filt(1:end-1,1:end-1),ts,rm);
% % % alpha_total = fillmissing(alpha_total,'nearest',1); % Fill in NaNs in dimension 1
% % % alpha_total = fillmissing(alpha_total,'nearest',2); % Fill in NaNs in dimension 2
% 
% T_final_test_pad = padarray(T_final_test,[8/2,30/2],'replicate');
% T_final_test_filt = filter2(k,T_final_test_pad,'valid');
% T_final_test_filt = interp2(ts-t_step/2,rm-rangeBin/2,T_final_test_filt(1:end-1,1:end-1),ts,rm);
% T_finalm = T_final_test_filt .* SNRm .* cloud_SDm_above;
% T_finalm(T_finalm <= 0) = NaN;

T_finalm = T_final_test .* SNRm .* cloud_SDm_above;
T_finalm(T_finalm <= 0) = NaN;

% % end %full loop


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




T_finalm = T_final_tests .* SNRm .* cloud_SDm_above ;
T_finalm(T_finalm <= 0) = NaN;


%%
disp('Gathering other data')

InsideCell = interp1(rmmissing(DataStructure2.Thermocouple.InsideCell.TimeStamp),rmmissing(DataStructure2.Thermocouple.InsideCell.Temperature),thr);
OutsideCell = interp1(rmmissing(DataStructure2.Thermocouple.InsideCell.TimeStamp),rmmissing(DataStructure2.Thermocouple.OutsideCell.Temperature),thr);
TSOA = interp1(rmmissing(DataStructure2.Thermocouple.InsideCell.TimeStamp),rmmissing(DataStructure2.Thermocouple.TSOA.Temperature),thr);
RoomTemp = interp1(rmmissing(DataStructure2.Thermocouple.InsideCell.TimeStamp),rmmissing(DataStructure2.Thermocouple.RoomTemp.Temperature),thr);

OnlineWaveDiff = interp1(rmmissing(DataStructure2.Laser.O2Online.TimeStamp  ),rmmissing(DataStructure2.Laser.O2Online.WaveDiff),thr);
OfflineWaveDiff = interp1(rmmissing(DataStructure2.Laser.O2Offline.TimeStamp  ),rmmissing(DataStructure2.Laser.O2Offline.WaveDiff),thr);

OnlineTemperatureActual = interp1(rmmissing(DataStructure2.Laser.O2Online.TimeStamp  ),rmmissing(DataStructure2.Laser.O2Online.TemperatureActual),thr);
OnlineTemperatureDesired = interp1(rmmissing(DataStructure2.Laser.O2Online.TimeStamp  ),rmmissing(DataStructure2.Laser.O2Online.TemperatureDesired),thr);
OfflineTemperatureActual = interp1(rmmissing(DataStructure2.Laser.O2Offline.TimeStamp  ),rmmissing(DataStructure2.Laser.O2Offline.TemperatureActual),thr);
OfflineTemperatureDesired = interp1(rmmissing(DataStructure2.Laser.O2Offline.TimeStamp  ),rmmissing(DataStructure2.Laser.O2Offline.TemperatureDesired),thr);

EtalonTemperatureActual = interp1(rmmissing(DataStructure2.Etalon.O2Etalon.TimeStamp  ),rmmissing(DataStructure2.Etalon.O2Etalon.TemperatureActual),thr);

channel0time = DataStructure2.MCS.Channel0.NewTimeGrid;
channel2time = DataStructure2.MCS.Channel2.NewTimeGrid;
channel8time = DataStructure2.MCS.Channel8.NewTimeGrid;
channel10time = DataStructure2.MCS.Channel10.NewTimeGrid;
channel0Bins = DataStructure2.MCS.Channel0.NBins;
channel2Bins = DataStructure2.MCS.Channel2.NBins;
channel8Bins = DataStructure2.MCS.Channel8.NBins;
channel10Bins = DataStructure2.MCS.Channel10.NBins;

%%
% savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\7_2_20\';
% %savePath = 'C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\7_2_20\';
% % savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\10_28_20\';
savePath = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\Reports\2_28_21\';
savePath = 'C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\2_28_21\';

savePath = 'C:\Users\d98b386\OneDrive - Montana State University - Bozeman\research s19\Reports\3_26_21\';
savetrue = 0;

if savetrue ==1
[y,m,d]=ymd(span_days(1));
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

% save([savePath date],'alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','o2on','o2off','o2on_mol','o2off_mol'...
%     ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
%     ,'absorption_f','InsideCell','OutsideCell','TSOA','RoomTemp','OnlineWaveDiff','OfflineWaveDiff'...
%     ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
%     ,'date_ts','sonde_ind','WV','Ts','exclusion','T_final_test','T_sonde','absorption_sonde')

save([savePath date 'nowv'],'alpha_totalm','alpha_0','T','T_finalm','BSR','ts','cloud_SDm','rm','o2on','o2off','o2on_mol','o2off_mol'...
    ,'Ts_fit','L_fit_sm_test','bg_o2off','bg_o2on','bg_o2off_mol','bg_o2on_mol','alpha_total_raw','absorption'...
    ,'InsideCell','OutsideCell','TSOA','RoomTemp','OnlineWaveDiff','OfflineWaveDiff'...
    ,'o2on_intp','o2off_intp','o2on_intp_mol','o2off_intp_mol'...
    ,'date_ts','sonde_ind','WV','Ts','exclusion','T_final_test','T_sonde','absorption_sonde');

% AerosolBackscatterCoefficient=LidarData.AerosolBackscatterCoefficient;
% MolecularBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient;
% save([savePath date],'BSR','rm','ts','Ts','Ps','MolecularBackscatterCoefficient','AerosolBackscatterCoefficient');
end
%%
% % % disp('calulating model photons')
% % % tic
% % % [N_on,N_off,N_onO, N_offO,N_off_mol,N_on_pulse,N_off_pulse,N_on_pulse_MCS,N_off_pulse_MCS,modelRange] = modelReturns(absorption_f,T_etalon_on,T_etalon_off,rm,ts,WV,lambda_online,lambda_offline,nu_scan_3D_short,...
% % %     online_index(1),offline_index(1),Ba,Bm,g1_m_on,g1_m_off,g1_a_on,g1_a_off);
% % % 
% % % ind_r_lo = 1:length(modelRange)-1;
% % % ind_r_hi = 2:length(modelRange);
% % % modelRangeBin = modelRange(2)-modelRange(1);
% % % 
% % % ln_o2 = log((N_on(ind_r_lo,:).*N_off(ind_r_hi,:))./(N_on(ind_r_hi,:).*N_off(ind_r_lo,:))); % Natural log of counts
% % % alpha_0_raw_model = ln_o2./2./(modelRangeBin);                              %[1/m] 
% % % %alpha_0_raw_model = ln_o2./2./(rangeBin);                              %[1/m] 
% % % ln_o2 = log((N_on_pulse(ind_r_lo,:).*N_off_pulse(ind_r_hi,:))./(N_on_pulse(ind_r_hi,:).*N_off_pulse(ind_r_lo,:))); % Natural log of counts
% % % alpha_0_raw_model_pulse = ln_o2./2./(modelRangeBin);                              %[1/m] 
% % % toc
%%

% alpha_0_pad = padarray(alpha_0_raw_model,[oversample/2,t_avg/2],'replicate');
% alpha_0_filt = filter2(k,alpha_0_pad,'valid');
% %alpha_0 = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2+rangeBin*oversample/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% alpha_0_model = interp2(ts-t_step/2,modelRange(ind_r_lo)-modelRangeBin/2,alpha_0_filt(1:end-1,1:end-1),ts,modelRange);
% %alpha_0 = interp2(ts,rm(ind_r_lo),alpha_0_raw,ts,rm);
% alpha_0_model = fillmissing(alpha_0_model,'nearest',1);                                 % Fill in NaNs in dimension 1
% alpha_0_model = fillmissing(alpha_0_model,'nearest',2);                                 % Fill in NaNs in dimension 2
% 
% alpha_0_pad = padarray(alpha_0_raw_model_pulse,[oversample/2,t_avg/2],'replicate');
% alpha_0_filt = filter2(k,alpha_0_pad,'valid');
% %alpha_0 = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2+rangeBin*oversample/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% alpha_0_model_pulse = interp2(ts-t_step/2,rm(ind_r_lo)-rangeBin/2,alpha_0_filt(1:end-1,1:end-1),ts,rm);
% %alpha_0 = interp2(ts,rm(ind_r_lo),alpha_0_raw,ts,rm);
% alpha_0_model_pulse = fillmissing(alpha_0_model_pulse,'nearest',1);                                 % Fill in NaNs in dimension 1
% alpha_0_model_pulse = fillmissing(alpha_0_model_pulse,'nearest',2);                                 % Fill in NaNs in dimension 2



% alpha_0_fit = zeros(i_range,i_time);
%  for i = 1:i_time
%      alpha_0_fit_obj = fit(rm(26:106,1),alpha_0_model(26:106,i),'poly1');
%      alpha_0_fit(:,i) = alpha_0_fit_obj(rm);
%  end

% alpha_1_raw = 0.5.*(alpha_0_raw_model.*W1 + G1);  %[1/m]
%  alpha_1 = alpha_1_raw;
%  
%  Tm1 = exp(-cumtrapz(rm,oversample.*alpha_1.*f,1));      %[none] First order transmission
% 
% % === Second Order ===
% % Integrated terms
% zeta_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1),3)*nuBin*100;             %[none]
% eta_Tm1_int = trapz(eta.*Tm0.*(1-Tm1),3)*nuBin*100;               %[1/m]
% zeta_ls_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1).*(1-f),3)*nuBin*100;   %[none]
% 
% W2 = (zeta_ls_int.*zeta_ls_Tm1_int./(zeta_int.^2)) - (zeta_ls_Tm1_int./zeta_int);   %[none]
% G2 = (eta_int.*zeta_Tm1_int./(zeta_int.^2)) - (eta_Tm1_int./zeta_int);              %[1/m]
% 
% alpha_2_raw = 0.5.*(alpha_1.*W1 + alpha_0_fit.*W2 + G2);
% 
% alpha_total_raw_model = alpha_0_model + alpha_1 + alpha_2;
 

%%

%===============
%=== Figures ===
%===============
% plot_time = 337;                         %[min];
% [~,p_point] = min(abs(plot_time-ts/60)); % Find closest value to 338min for comparison to other program
                        %[hr];
plot_time=1.0833;
%plot_time=1.5833;
%plot_time=2.5833;
%plot_time=10;
[~,p_point] = min(abs(plot_time-ts/60/60)); % Find closest value to 338min for comparison to other program
p_point(1:length(rm),1)=p_point;

sonde_index =1;
p_point = sonde_ind(:,sonde_index);

figure(728590)
plot(ts/60/60,permute(L_fit_sm_test(1,:,end)*1000,[3 2 1]))
hold on
plot(ts([p_point(1) p_point(end)])/60/60,[-10 -4])
hold off
yline(-6.5);
grid on
title('Fitted lapse rate')
title(sprintf('Fitted lapse rate\n %s',span_days(1)))
xlabel('Time (UTC hr)')
ylabel('Lapse rate (K/km)')

figure(728591)
%plot(ts/60/60,Ts - permute(Ts_fit(1,:,:),[3 2 1]))
plot(ts/60/60,T(1,:) - permute(Ts_fit(1,:,:),[3 2 1]))
hold on
plot(ts([p_point(1) p_point(end)])/60/60,[-20 20])
hold off
grid on
title('Surface temp - fitted surface temp')
xlabel('Time (UTC hr)')
ylabel('\deltaT (K)')

figure(9985)
semilogy(thr,bg_o2off)
hold on
semilogy(thr,bg_o2on)
semilogy(thr,bg_o2off_mol)
semilogy(thr,bg_o2on_mol)
title('Background')
xlabel('Time UTC hr')
ylabel('photons')

plot(ts([p_point(1) p_point(end)])/60/60,[1 100])
hold off
legend('Offline','Online','Offline molecular','Online molecular')
grid on


% figure(9986)
% semilogy(thr,bg_o2off/(nsPerBin*summedBins)*10^-6)
% hold on
% semilogy(thr,bg_o2on/(nsPerBin * summedBins)*10^-6)
% semilogy(thr,bg_o2off_mol/(nsPerBin * summedBins)*10^-6)
% semilogy(thr,bg_o2on_mol/(nsPerBin * summedBins)*10^-6)
% title('Background count rate')
% xlabel('Time UTC hr')
% ylabel('photoelectrons/\mu s')
% xline(thr(p_point));
% hold off
% legend('Offline','Online','Offline molecular','Online molecular')
% grid on

% figure(9986)
% semilogy(thr,bg_o2off)
% hold on
% semilogy(thr,bg_o2on)
% semilogy(thr,bg_o2off_mol)
% semilogy(thr,bg_o2on_mol)
% title('Background counts')
% xlabel('Time UTC hr')
% ylabel('photoelectrons/\mu s')
% xline(thr(p_point));
% hold off
% legend('Offline','Online','Offline molecular','Online molecular')
% grid on

% % figure(118)
% % [~,online_index1] = min(abs( lambda_scan_3D_short(1,1,:) - 769.7958));
% % [~,online_index2] = min(abs( lambda_scan_3D_short(1,1,:) - 769.7948));
% % [~,online_index3] = min(abs( lambda_scan_3D_short(1,1,:) - 769.7968));
% % %plot(permute(del_nu(1,1,:),[3 2 1]),permute(absorption_f(:,p_point,:),[3 1 2]))
% % xline(769.7958);
% % xline(769.7968);
% % xline(769.7948);
% % hold on
% % plot(permute(lambda_scan_3D_short(1,1,:),[3 2 1]),permute(absorption_f(:,p_point,:),[3 1 2]))
% % hold on
% % %plot(permute(del_nu(1,1,:),[3 2 1]),permute(absorption_f2(:,p_point,:),[3 1 2]))
% % xlabel('Wavelength (nm)')
% % ylabel('Absorption (m^{-1})')
% % legend('769.7958','769.7968','769.7948')
% % 
% % %xline(permute(del_nu(1,1,online_index1(1)),[3 2 1]));
% % %xline(permute(del_nu(1,1,online_index2(1)),[3 2 1]));
% % grid on
% % hold off
% % %legend('absorption_f','absorption_f2')
% % title('O2 absorption spectrum')




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

% % figure(99588)
% % semilogy(permute(lambda_scan_3D_short(:,:,:),[3 2 1]),permute(doppler_O2_ret(:,p_point,:)./max(doppler_O2_ret(:,p_point,:),[],3),[3 1 2]))
% % hold on
% % semilogy(permute(lambda_scan_3D_short(:,:,:),[3 2 1]),permute(T_etalon_on(:,:,:),[3 1 2]))
% % 
% % semilogy(permute(lambda_scan_3D_short_off(:,:,:),[3 2 1]),permute(doppler_O2_ret_off(:,p_point,:)./max(doppler_O2_ret_off(:,p_point,:),[],3),[3 1 2]))
% % semilogy(permute(lambda_scan_3D_short_off(:,:,:),[3 2 1]),permute(T_etalon_off(:,:,:),[3 1 2]))
% % xline(lambda_offline(1));
% % xline(lambda_online(1));
% % hold off
% % grid on
% % ylim([10^-4 1])
% % 
% % xlabel('Wavelength')
% % title('doppler broadending')

% figure(99589)
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(g1(:,p_point,:),[3 1 2]),'-*')
% grid on
% xlabel('wavenumber')
% title('backscatter lineshape')

figure(99492)
plot(diag(BSR(:,p_point)),rkm)
%title('Bacscatter ratio (Bt/Bm)')
%title(sprintf('Bacscatter ratio (Bt/Bm)\n %s at %.2f UTC',span_days(1),thr(p_point)))
if (thr(p_point(1))-floor(thr(p_point(1))))*60 < 10
    title(sprintf('Bacscatter ratio (Bt/Bm)\n %s at %.0f:0%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
else
    title(sprintf('Bacscatter ratio (Bt/Bm)\n %s at %.0f:%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
end
grid on
ylabel('Range (km)')

figure(99493)
plot(diff(BSR(:,p_point))./diff(rm),rkm(1:end-1)+(rkm(2)-rkm(1))/2)
title('Bacscatter ratio d(Bt/Bm)/dr')
grid on
ylabel('Range (km)')


% figure(99493)
% 
%     [~,g] = sgolay(2,oversample);
%     parfor j=1:i_time
%         dBSRSG(:,j) = conv(BSR(:,j), factorial(1)/(-rangeBin)^1 * g(:,2), 'same');
%     end
% plot(diag(dBSRSG(:,p_point)),rkm)
% hold on
% dBSR = (BSR(ind_r_hi,:)-BSR(ind_r_lo,:))/(rangeBin*oversample);
% %plot(diag(dBSR(:,p_point)),rkm(ind_r_lo+4))
% plot(diag(dBSR(:,p_point)),rkm(ind_r_lo))
% 
% dBSRlo = diff(BSR)/rangeBin;
% plot(diag(dBSRlo(:,p_point)),rkm(1:end-1))
% hold off
% legend('2SG','over forward','forward')
% %title({'dBSR/dr';datestr(date_ts_N(p_point(1)))})
% grid on
% ylabel('Range (km)')

% % figure(99494)
% % plot(diag(alpha_0(:,p_point)-absorption_sonde{sonde_index}),rkm)
% % 
% % %close(99495)
% % figure(99495)
% % cla reset
% % line(diag(dBSRSG(:,p_point)),rkm,'Color','r')
% % line(diag(dBSRSG(:,p_point)./BSR(:,p_point)),rkm,'Color','b')
% % %line(diag(dBSRlo(:,p_point)),rkm(1:end-1))
% % legend('dBSR/dr')
% % xlim([-.01 .01])
% % xlabel('dBSR/dr')
% % ylabel('Range (km)')
% % 
% % ax1 = gca; % current axes
% % ax1.XColor = 'r';
% % ax1.YColor = 'r';
% % 
% % ax1_pos = ax1.Position; % position of first axes
% % ax2 = axes('Position',ax1_pos,...
% %     'XAxisLocation','top',...
% %     'YAxisLocation','right',...
% %     'Color','none');
% % 
% % line(diag(alpha_0(:,p_point)-absorption_sonde{sonde_index}),rkm,'Parent',ax2,'Color','k')
% % 
% % grid on
% % xlim([-.5 .5]*10^-3)
% % xlabel('\alpha_0 - \alpha_{sonde}')
% % ylabel('Range (km)')
% % legend('\alpha_0 - \alpha_{sonde}')

% figure(99595)
% %plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:),[3 1 2]),'-*')
% hold on
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:).*(1-f(:,p_point,:)),[3 1 2]),'-*')
% hold off
% grid on
% xlabel('wavenumber')
% title('zeta1s')
% 
% figure(99596)
% %plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:).*Tm0(:,p_point,:),[3 1 2]),'-*')
% hold on
% plot(permute(nu_scan_3D_short,[3 2 1]),permute(zeta(:,p_point,:),[3 1 2]),'-*')
% hold off
% grid on
% xlabel('wavenumber')
% title('zeta1')

figure(99499)
subplot(2,1,1)
plot(diag(WV(:,p_point)),rm)
xlabel('WV (molec/m^3')
ylabel('Range (m)')
legend('VW')
subplot(2,1,2)
%plot(sondeStruc.AH,(sondeStruc.Height-sondeStruc.Height(1))/1000)
xlabel('Absolute Humidity g/m^3')

% figure(994990)
% plot(WV(:,p_point)*18.01528*6.022e-23,rm)
% xlabel('WV (g/m^3')
% ylabel('Range (m)')
% legend('VW')



figure(884)
%plot(diag(T(:,p_point)),rkm)
%plot(T(:,p_point)+2,rkm)
%plot(T(:,p_point)-2,rkm)
%plot(Ts(p_point),0,'+')
%plot(Ts_fit(1,p_point,end),0,'+')

plot(diag(T_finalm(:,p_point)),rkm)
hold on
plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end)),rkm,'--')
%plot(T_final(:,p_point),rm)
plot(diag(T(:,p_point)),rkm)
exclusion(exclusion==0)=NaN;
plot(diag(T_final_test(:,p_point)).*diag(exclusion(:,p_point,end)),rkm,'*')
plot(T_sonde(:,sonde_index),rm/1000,'.-')
%plot(Tnoaa(:,p_point),rmnoaa(:,p_point)/1000)
hold off
grid on
xlabel('Temperature (K)')
ylabel('Range (km)')
ylim([0 6])
%title('temp guess and retrieved from pertabative absorption')
%title(sprintf('Temp guess and retrieved from pertabative absorption\n %s at %.2f UTC',span_days(1),thr(p_point)))
if (thr(p_point(1))-floor(thr(p_point(1))))*60 < 10
    title(sprintf('Temp guess and retrieved from pertabative absorption\n %s at %.0f:0%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
else
    title(sprintf('Temp guess and retrieved from pertabative absorption\n %s at %.0f:%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
end
%legend('Temp guess','Retrieved Temperature','Fitted lapse rate','Points excluded from fit','Sonde temperature','Location','southwest')
%legend('Retrieved Temperature','Sonde temperature','Location','southwest')
legend('Retrieved Temperature','Fitted temp','temp guess','exclusion','Sonde temperature','Location','southwest')
%legend('Retrieved Temperature','Fitted temp','temp guess','exclusion','Location','southwest')

figure(885)
plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end))-diag(T_finalm(:,p_point)),rkm)

%title(sprintf('Temperature Deviation  %s (UTC)',date_time))
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-1 -1],[0 6],'Color','k','LineStyle','-.')
line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color','k')
line([2 2],[0 6],'Color','k')
xlabel('\Delta T (T_{fit} - T_{retrieved}) (K)')
ylabel('Range (km)')
legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','northwest')
%legend(sprintf('%s at %.2f UTC',span_days(1),thr(p_point(1))),sprintf('%s at %.2f UTC',span_days(1),thr(sonde_ind(1,2))))
ylim([0 6])
xlim([-10 10])
%title(sprintf('Temperature difference\n %s at %.2f UTC',span_days(1),thr(p_point(1))))
hold off

figure(888)
%plot(T(:,p_point) - T_finalm(:,p_point),rkm)
plot(T_sonde(:,sonde_index) - diag(T_finalm(:,p_point)),rkm)
hold on
%%%%%plot(T_sonde(:,2)-diag(T_finalm(:,sonde_ind(:,2))),rkm,'--')
%title(sprintf('Temperature Deviation  %s (UTC)',date_time))
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-1 -1],[0 6],'Color','k','LineStyle','-.')
line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color','k')
line([2 2],[0 6],'Color','k')
%xlabel('\Delta T (T_{model} - T_{retrieved}) (K)')
hold off
xlabel('\Delta T (T_{sonde} - T_{retrieved}) (K)')
ylabel('Range (km)')
%%%%%%%%%%%%legend(sprintf('%s at %.2f UTC',span_days(1),thr(p_point(1))),sprintf('%s at %.2f UTC',span_days(1),thr(sonde_ind(1,2))))
%legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','Southwest')
ylim([0 5])
xlim([-10 10])
%title(sprintf('Temperature difference\n %s at %.2f UTC',span_days(1),thr(p_point)))
if (thr(p_point(1))-floor(thr(p_point(1))))*60 < 10
    title(sprintf('Temperature difference\n %s at %.0f:0%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
else
    title(sprintf('Temperature difference\n %s at %.0f:%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
end




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
[~,ind_km_max] = min(abs(rkm-r_max_plot)); % Max range Index
imAlpha=ones(size(T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
colorLimits = [min(T_finalm(1:ind_km_max,:)-1,[],'all'), max(T_finalm(1:ind_km_max,:),[],'all')]; % Set color limits to only include data -1
imagesc(thr,rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits) %Plot
hold on
%sonde line
plot(ts([p_point(1) p_point(end)])/60/60,[rkm(1) rkm(end)])
colors = colormap; % Set colors to colormap
colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colormap(colors); % Make colormap new colors
set(gca, 'YDir','normal') % Set y axis increasing
set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
title(sprintf('Temperature retrieval\n%s(UTC)',span_days(1))) %Add title
xlabel('Time (UTC hours)')
ylabel('Range (km)')
title(colorbar,'Temperature (K)') % Add title to colorbar
hold off

%xline(thr(sonde_ind(2)))


% % % % % % figure(887)
% % % % % % r_max_plot = 6; % Max range to plot [km]
% % % % % % [~,ind_km_max] = min(abs(rkm-r_max_plot));
% % % % % % imAlpha=ones(size(alpha_totalm(1:ind_km_max,:)));
% % % % % % imAlpha(isnan(alpha_totalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
% % % % % % imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
% % % % % % %imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1;
% % % % % % colorLimits = [min(alpha_totalm(1:ind_km_max,:),[],'all'), max(alpha_totalm(1:ind_km_max,:),[],'all')];
% % % % % % imagesc(thr,rkm(1:ind_km_max),alpha_totalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% % % % % % hold on
% % % % % % xline(thr(p_point),'color','b');
% % % % % % colors = colormap;
% % % % % % colors(1,:) = [0 0 0];%set lowest color black
% % % % % % colormap(colors);
% % % % % % set(gca, 'YDir','normal')
% % % % % % set(gca,'color',[1 1 1]);%color background white
% % % % % % colorbar
% % % % % % %title(sprintf('Temperature %s to %s (UTC)',ts_datetime(1),ts_datetime(end)))
% % % % % % title(sprintf('Alpha final retrieval on %s(UTC)',span_days(1)))
% % % % % % xlabel('Time (UTC hours)')
% % % % % % ylabel('Range (km)')
% % % % % % title(colorbar,'Absorption (m^-1)')
% % % % % % hold off






figure(4)
plot(diag(alpha_0(:,p_point)),rkm)
hold on
%%plot(absorption(:,p_point),rkm,'--')

% % [~,online_index1] = min(abs( lambda_scan_3D_short(1,1,:) - 769.7968));
% % [~,online_index2] = min(abs( lambda_scan_3D_short(1,1,:) - 769.7948));
% % plot(permute(absorption_f(:,p_point,online_index1:online_index2),[3 1 2]),rkm)
%plot(rm,absorption_new(:,online_index_new)-absorption_new(:,offline_index_new))

plot(diag(alpha_totalm(:,p_point)),rkm,'LineWidth',2)
%plot(alpha_0m(:,p_point)+alpha_1(:,p_point),rkm)
%%%
%plot(alpha_total(:,p_point),rkm)
%%%plot(alpha_0_raw_model(:,p_point),modelRange(ind_r_lo)/1000,'--')
%%%plot(alpha_0_raw_model_pulse(:,p_point),modelRange(ind_r_lo)/1000,'--')
plot(diag(alpha_total_raw(:,p_point)),rkm)
plot(absorption_sonde{sonde_index},rkm,'.-')
%plot(alpha_total_raw_model(:,p_point),rkm,'--')
%plot(alpha_0_fit(:,p_point),rkm)
%plot(alpha_0(:,p_point)+d_alpha(:,p_point),rkm,'.-')

%plot(diag(alpha_0_mol(:,p_point)),rkm,'--')
%plot(diag(alpha_0_mol_corr(:,p_point)),rkm,'.-')

%plot(diag(alpha_total(:,p_point)),rkm,'.-')
hold off
grid on
%title('Calculated absorption coefficient')
if (thr(p_point(1))-floor(thr(p_point(1))))*60 < 10
    title(sprintf('Calculated absorption coefficient\n %s at %.0f:0%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
else
    title(sprintf('Calculated absorption coefficient\n %s at %.0f:%.0f UTC',span_days(1),floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60))
end
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')
xlim([-0.5e-4 4e-4])
%%%%legend('Measured zeroth order absorption','Theoretical absorption online','Pertabative absorption','Alpha total raw','Alpha raw model','alpha raw pulse model')
%legend('Measured zeroth order absorption','Theoretical absorption online','Pertabative absorption','Alpha total raw','Sonde Absorption')
%legend('Measured zeroth order absorption','Pertabative absorption','Sonde Absorption')
%legend('Measured zeroth order absorption','Pertabative absorption','Alpha 0+1+2','Sonde')%,'theoretical absorption online-offine')
legend('Measured zeroth order absorption','Pertabative absorption','Alpha 0+1+2','Radiosonde model')


figure(7)
plot(diag(o2on(:,p_point)),rm)
hold on
plot(diag(o2off(:,p_point)),rm)
plot(diag(o2on_mol(:,p_point)),rm)
plot(diag(o2off_mol(:,p_point)),rm)

% % plot(modelRange,N_on(:,p_point),'.-')
% % plot(modelRange,N_off(:,p_point),'.-')
% % plot(modelRange,N_off_mol(:,p_point),'.-')

plot(diag(o2on_noise(1:i_range,p_point)),rm,'--')
plot(diag(o2off_noise(1:i_range,p_point)),rm,'--')
plot(diag(o2on_noise_mol(1:i_range,p_point)),rm,'--')
plot(diag(o2off_noise_mol(1:i_range,p_point)),rm,'--')

%%%%%%plot(diag(o2off_mol_corr(:,p_point)),rm,'.-')


xlim([0 500])
%view([90 -90])
title('Averaged photon returns')
legend('Online','Offline','Online molecular','Offline molecular','raw online','raw offline','raw online molecular','raw offline molecular')
hold off
grid on
xlabel('Range [m]')
ylabel('Photons')

% % figure(8)
% % plot(rm,N_on(:,p_point)-o2on(:,p_point))
% % hold on
% % plot(rm,N_off(:,p_point)-o2off(:,p_point))
% % plot(rm,N_off_mol(:,p_point)-o2off_mol(:,p_point))
% % ylim([-30 30])
% % legend('online','offline','offline molecular')
% % title('Model photon returns - averaged photon counts')
% % xlabel('Range (m)')
% % ylabel('Model - collected (counts)')
% % hold off


% figure(8)
% plot(rm,o2on(:,p_point).*rm.^2)
% hold on
% plot(rm,o2off(:,p_point).*rm.^2)
% plot(rm,o2on_mol(:,p_point).*rm.^2)
% plot(rm,o2off_mol(:,p_point).*rm.^2)
% 
% plot(rm,o2on(:,p_point)./o2on_correctionFactor(1:i_range,p_point).*rm.^2)
% hold on
% plot(rm,o2off(:,p_point)./o2off_correctionFactor(1:i_range,p_point).*rm.^2)
% plot(rm,o2on_mol(:,p_point)./o2on_correctionFactor_mol(1:i_range,p_point).*rm.^2)
% plot(rm,o2off_mol(:,p_point)./o2off_correctionFactor_mol(1:i_range,p_point).*rm.^2)
% title('Averaged photon returns multiplied by r^2')
% legend('Online','Offline','Online molecular','Offline molecular')
% xlabel('Range m')
% ylabel('Photons * r^2')
% grid on
% hold off

% figure(9)
% plot(rm,o2on_correctionFactor(1:i_range,p_point))
% hold on
% plot(rm,o2off_correctionFactor(1:i_range,p_point))
% plot(rm,o2on_correctionFactor_mol(1:i_range,p_point))
% plot(rm,o2off_correctionFactor_mol(1:i_range,p_point))
% ylim([-inf 1.1])
% title('ApdCorrection factor')
% legend('Online','Offline','Online molecular','Offline molecular')
% xlabel('Range m')
% ylabel('Unitless')
% grid on
% hold off

figure(738732)
subplot(4,1,1)
plot(thr,InsideCell)
legend('Inside Cell')
title('Thermocouples')
ylabel('^o C')
subplot(4,1,2)
plot(thr,OutsideCell)
ylabel('^o C')
legend('Outside Cell')
subplot(4,1,3)
plot(thr,TSOA)
ylabel('^o C')
legend('TSOA')
subplot(4,1,4)
plot(thr,RoomTemp)
ylabel('^o C')
legend('Room Temp')

xlabel('Time UTC hr')



% figure(738733)
% plot(thr,DataStructure2.Laser.O2Online.TemperatureActual)
% hold on
% plot(thr,DataStructure2.Laser.O2Offline.TemperatureActual)
% hold off
% legend('Online','Offline')
% title('wavelength Difference')

figure(738734)
plot(thr,OnlineWaveDiff)
hold on
plot(thr,OfflineWaveDiff)
hold off
legend('Online','Offline')
title('wavelength Difference')
ylim([-5*10^-4 5*10^-4])
ylabel('(nm)')
xlabel('Time hr')

figure(738735)
subplot(3,1,1)
plot(thr,OnlineTemperatureActual)
hold on
plot(thr,OnlineTemperatureDesired)
legend('Actual','Desired')
hold off
title('Online Laser Temperature')
ylabel('(C)')
subplot(3,1,2)
plot(thr,OfflineTemperatureActual)
hold on
plot(thr,OfflineTemperatureDesired)
legend('Actual','Desired')
hold off
title('Offline Laser Temperature')
ylabel('(C)')
subplot(3,1,3)
plot(thr,EtalonTemperatureActual)
title('Etalon Temperature')
ylabel('(C)')
xlabel('Time hr')


figure(9999)
plot(rm_raw_o2/1000,o2on_intp(:,p_point))
hold on
plot(rm_raw_o2/1000,o2off_intp(:,p_point))
plot(rm_raw_o2/1000,o2off_intp_mol(:,p_point))
plot(rm_raw_o2/1000,o2on_intp_mol(:,p_point))
hold off
title('Raw counts')
ylabel('Counts')
legend('On Combined','Off Combined','On Molecular','Off Molecular')
xlabel('Range (km)')
ylim([0 500])
xlim([-1 15])
grid on


% % figure(999999)
% % plot_time=23.6667;
% % [~,p_point1] = min(abs(plot_time-ts/60/60)); % Find closest value to 338min for comparison to other program
% % plot_time=23.75;
% % [~,p_point2] = min(abs(plot_time-ts/60/60)); % Find closest value to 338min for comparison to other program
% % pulsePoint=8;
% % 
% % 
% % apdx=rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000;
% % apdy=mean(o2on_bgsub(pulsePoint:end,p_point1+2:p_point2+1),2);
% % [~,p_point3] = min(abs(5000-rm_raw_o2)); % Find closest value to 338min for comparison to other program
% % apdx=rm_raw_o2(pulsePoint:p_point3)/1000-rm_raw_o2(pulsePoint)/1000;
% % apdy=mean(o2on_intp(pulsePoint:p_point3,p_point1+2:p_point2+1),2);
% % apdy2=mean(o2off_intp(pulsePoint:p_point3,p_point1+2:p_point2+1),2);
% % apdy3=mean(o2on_intp_mol(pulsePoint:p_point3,p_point1+2:p_point2+1),2);
% % 
% % %plot(rm_raw_o2/1000,o2on_intp(:,p_point1:p_point2))
% % plot(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000,mean(o2on_intp(pulsePoint:end,p_point1+2:p_point2+1),2))
% % %o2onFit = fit(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000,mean(o2on_bgsub(pulsePoint:end,p_point1+2:p_point2+1),2),'smoothingspline');
% % 
% % hold on
% % plot(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000,18.84*exp(-(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000)/0.1382)+3.838)
% % 
% % %plot(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000,mean(o2off_intp(pulsePoint:end,p_point1+2:p_point2+1),2))
% % %plot(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000,mean(o2on_intp_mol(pulsePoint:end,p_point1+2:p_point2+1),2))
% % %plot(rm_raw_o2(pulsePoint:end)/1000-rm_raw_o2(pulsePoint)/1000,mean(o2off_intp_mol(pulsePoint:end,p_point1+2:p_point2+1),2))
% % 
% % %plot(rm_raw_o2/1000,o2off_intp(:,p_point1:p_point2))
% % %plot(rm_raw_o2/1000,o2off_intp_mol(:,p_point1:p_point2))
% % %plot(rm_raw_o2/1000,o2on_intp_mol(:,p_point1:p_point2))
% % hold off
% % title('Raw counts')
% % ylabel('Counts')
% % legend('On Combined','Off Combined','On Molecular','Off Molecular')
% % xlabel('Range (km)')
% % ylim([0 50])
% % xlim([-1 15])
% % grid on



% figure(99999)
% plot(rm,o2on_noise(:,p_point))
% hold on
% plot(rm,o2off_noise(:,p_point))
% plot(rm,o2off_noise_mol(:,p_point))
% plot(rm,o2on_noise_mol(:,p_point))
% 
% %plot(rm_raw_o2,A_co*exp(-rm_raw_o2/1000/RC_co))
% %plot(rm_raw_o2,o2on_bgsub(1,p_point).*exp(-rm_raw_o2/1000/RC_co))
% 
% plot(rm,smooth(o2on_noise(:,p_point),10,'loess'),'--')
% plot(rm,smooth(o2off_noise(:,p_point),10,'loess'),'--')
% plot(rm,smooth(o2off_noise_mol(:,p_point),10,'loess'),'--')
% plot(rm,smooth(o2on_noise_mol(:,p_point),10,'loess'),'--')
% hold off
% title('raw counts')

% figure(783133)
% changeWavelength=0.0005;%[nm]
% [~,onlineplus] = min(abs(lambda_scan(online_index(p_point))+changeWavelength-lambda_scan));
% [~,onlineminus] = min(abs(lambda_scan(online_index(p_point))-changeWavelength-lambda_scan));
% 
% plot(lambda_scan,permute(absorption_f(:,p_point,:),[1 3 2]))
% xline(lambda_scan(online_index(p_point)));
% xline(lambda_scan(onlineplus));
% xline(lambda_scan(onlineminus));
% xlabel('Wavelength')
% ylabel('Absorption')
% grid on

% % o2on_noise_sub = o2on_noise(:,p_point);
% o2on_noise_sub=smooth(o2on_noise(:,p_point),10,'loess');
% % o2off_noise_mol_sub = o2off_noise_mol(:,p_point);
% o2off_noise_sub=smooth(o2off_noise(:,p_point),10,'loess');
% % o2off_noise_sub = o2off_noise(:,p_point);
% o2off_noise_mol_sub=smooth(o2off_noise_mol(:,p_point),10,'loess');
% % o2on_noise_mol_sub = o2on_noise_mol(:,p_point);
% o2on_noise_mol_sub=smooth(o2on_noise_mol(:,p_point),10,'loess');
% save('noiseSub.mat','o2on_noise_sub','o2off_noise_mol_sub','o2off_noise_sub','o2on_noise_mol_sub')

% figure(738735)
% plot(thr,lambda_online)
% hold on
% plot(thr,lambda_offline)
% hold off
% legend('Online','Offline')
% title('wavelength Difference')

% figure(1999)
% plot(rm/1000,overlapFile)
% ylim([0 1.1])
% yline(1)
% legend('before','after','Location','southeast')
% title('Norm com online / norm mol online')

figure(88943)
plot(rm_raw_o2/1000,log(o2on_bgsub(:,p_point(1)).*rm_raw_o2.^2));
hold on
plot(rm_raw_o2/1000,log(o2off_bgsub(:,p_point(1)).*rm_raw_o2.^2));
plot(rm_raw_o2/1000,log(o2on_bgsub_mol(:,p_point(1)).*rm_raw_o2.^2));
plot(rm_raw_o2/1000,log(o2off_bgsub_mol(:,p_point(1)).*rm_raw_o2.^2));

%plot(rm/1000,log(o2on_noise(:,p_point(1)).*rm.^2),'--');
legend('on','off','on mol','off mol')
grid on
hold off



% % figure(889433)
% % plot(modelRange,log(N_on(:,p_point).*modelRange.^2))
% % hold on
% % plot(modelRange,log(N_off(:,p_point).*modelRange.^2))
% % plot(rm,log(o2on(:,p_point).*rm.^2),'--')
% % plot(rm,log(o2off(:,p_point).*rm.^2),'--')
% % plot(rm,log(o2on_mol(:,p_point).*rm.^2),'--')
% % plot(rm,log(o2off_mol(:,p_point).*rm.^2),'--')
% % 
% % legend('model on','model off','on','off','on mol','off mol')
% % hold off


figure(634432)
subplot(4,1,1)
bar(channel0time,channel0Bins)
title('Channel0')
subplot(4,1,2)
bar(channel2time,channel2Bins)
title('Channel2')
subplot(4,1,3)
bar(channel8time,channel8Bins)
title('Channel8')
subplot(4,1,4)
bar(channel10time,channel10Bins)
title('Channel10')

figure(83)
subplot(2,1,1)
plot(thr,Ts)
title('Weather Station Surface Temperature')
xlabel('Time')
ylabel('K')
subplot(2,1,2)
plot(thr,Ps)
title('Weather Station Suface Pressure')
xlabel('Time')
ylabel('Pressure ATM')


% figure()
% plot(o2off_mol(:,p_point)./(N_off_mol(1:2:end-1,p_point)*1),rm)
%close(648)
figure(648)
line(DataStructure2.UPS.all.TimeStamp,DataStructure2.UPS.all.BatteryTimeLeft,'Color','r')
line(DataStructure2.UPS.all.TimeStamp,DataStructure2.UPS.all.HoursOnBattery,'Color','g')
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('Battery remaining (hrs)')
legend('Battery remaining','Hours on battery','Location','Southwest')


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

line(DataStructure2.UPS.all.TimeStamp,DataStructure2.UPS.all.UPSTemperature,'Parent',ax2,'Color','k')
%xlabel('time (hrs)')
ylabel('Battery temperature (C)')



figure(589)
subplot(2,1,1)
imagesc(thr,rm/1000,o2on_noise.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
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
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
title('o2off unaveraged')

figure(5899)
subplot(2,2,1)
imagesc(thr,rm_raw_o2/1000,o2on_intp.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')
title('o2 on unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(thr,rm_raw_o2/1000,o2off_intp.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
title('o2off unaveraged')

subplot(2,2,3)
imagesc(thr,rm_raw_o2/1000,o2on_bgsub.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')
title('o2 on unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,4)
imagesc(thr,rm_raw_o2/1000,o2off_bgsub.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
title('o2off unaveraged bgsub')


figure(6899)
subplot(2,2,1)
imagesc(thr,rm_raw_o2/1000,o2on_intp_mol.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')
title('o2 on mol unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(thr,rm_raw_o2/1000,o2off_intp_mol.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
title('o2off mol unaveraged')

subplot(2,2,3)
imagesc(thr,rm_raw_o2/1000,o2on_bgsub_mol.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
%colormap('turbo')
title('o2 on mol unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,4)
imagesc(thr,rm_raw_o2/1000,o2off_bgsub_mol.*rm_raw_o2.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
%set(gca,'ColorScale','log')
colorbar
title('o2off mol unaveraged bgsub')


%Tick spacing
tickHours = 8;
tickMin = 0;
date2=datetime(2019,4,20,tickHours,tickMin,0);
date1=datetime(2019,4,20,0,0,0);
tickSpacing = datenum(date2)-datenum(date1);
tickAngle = 0;

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
% for i = 1:size(data_col_real,2)
%     if ~isnan(data_col_real(1,i))
%         %disp('ran')
%         %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
%         plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
%     end
% end
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
% for i = 1:size(data_col_real,2)
%     if ~isnan(data_col_real(1,i))
%         %disp('ran')
%         %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
%         plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
%     end
% end
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
% for i = 1:size(data_col_real,2)
%     if ~isnan(data_col_real(1,i))
%         %disp('ran')
%         %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
%         plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
%     end
% end
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
% for i = 1:size(data_col_real,2)
%     if ~isnan(data_col_real(1,i))
%         %disp('ran')
%         %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
%         plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
%     end
% end
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
% for i = 1:size(data_col_real,2)
%     if ~isnan(data_col_real(1,i))
%         %disp('ran')
%         %plot(ts([datenum(date_ts(data_col_real(1,i))) datenum(date_ts(data_col_real(end,i)))])/60/60,[rkm(1) rkm(end)])
%         plot(datenum(date_ts_N([data_col_real(1,i) data_col_real(end,i)])),[rkm(1) rkm(end)],'--k')
%     end
% end
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


BSRm = BSR .* SNRm .* cloud_SDm ;
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



