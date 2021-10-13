function [weather_Temperature_interp, weather_absPressure_interp, weather_WV_interp] = wunderWeather(span_days,ts,path)
%File: wunderWeather.m
%Date: 1/8/2021
%Author: Owen Cruikshank
%Inputs:
%   -span_days: datetime vector in UTC time of days. ex.: datetime(2020,2,22,'TimeZone','UTC');%yyyy,mm,dd
%   -ts:[s] main time grid of analysis
%   -path: string of file path to Weather station data folder
%
%Outputs:
%   -weather_Temperature_interp: [C] surface temperature from weather
%   station
%   -weather_absPressure_interp: [mbar] surface absolute pressure from weather station
%   -weather_VW_interp: [1/cm^3] surface water vapor number density

%Create datetime vector corresponding to time grid
thr = ts/60/60;%time in hours
time_hours = floor(thr);
time_minutes = floor((thr-time_hours).*60);
time_seconds = (((thr-time_hours).*60)-time_minutes)*60;
for i = 1:length(span_days(1,:))
    time_date(thr<24*i & thr>=24*i-24) = span_days(:,i);
    time_hours(thr<24*i & thr>=24*i-24) = time_hours(thr<24*i & thr>=24*i-24) - 24*(i-1);
end
time_date_hms = datetime(year(time_date),month(time_date),day(time_date),time_hours,time_minutes,time_seconds,'TimeZone','UTC');
% 
% %Seperate datetime vector to month and year
% month_span = month(span_days);
% if month_span < 10
%     month_span = cat(2,'0',num2str(month_span));
% else
%     month_span = num2str(month_span);
% end
% year_span = year(span_days);
% year_span = num2str(year_span);


%Read data in from file using month and year
file_date = 'wunderWeatherNovember.csv';
weatherFilePath = fullfile(path,'Weather station data',file_date);

% Read data
weatherFile = fileread(weatherFilePath);
% % Remove unwanted characters
% weatherFile = strrep(weatherFile,'''','');
% weatherFile = strrep(weatherFile,' K','');
% weatherFile = strrep(weatherFile,'"','');

% Read file into cell arrays
formatSpec='%s %s %s %s %s %s %s %s %s %s %s';
weather = textscan(weatherFile,formatSpec,'Delimiter',',','HeaderLines',1,'TreatAsEmpty','');

%create datetime vector for date from weather
weather_dateTime = [weather{1} weather{2}];

for i = 1:length(weather_dateTime(:,1))
    weather_dateTime2{i} = strjoin(weather_dateTime(i,:)');   
end
weather_dateTime = datetime(weather_dateTime2','TimeZone','America/Denver','InputFormat','M/d/yyyy h:mm a','Format','d-MMM-y HH:mm:ss Z');
%weather_dateTime = datetime(weather_dateTime2','TimeZone','UTC','InputFormat','M/d/yyyy h:mm a','Format','d-MMM-y HH:mm:ss Z');
%convert Denver to UTC
weather_dateTime.TimeZone = 'UTC';

%find data we want and convert to matricies
day_ = day(weather_dateTime) == day(span_days);
day_ = logical(sum(day_,2));
weather_dateTime = weather_dateTime(day_);

%weather_absPressure = cell2mat(weather(5)) .* 33.86389;%surface pressure [mbar]
for i = 1:length(weather{1,9})
    weatherPressure(i) = str2double(erase(weather{1,9}{i},[char(160) 'in'])).* 33.86389;
    weatherTemperature(i) = str2double(erase(weather{1,3}{i},[char(160) 'F']));
    weatherWV(i) = str2double(erase(weather{1,5}{i},[char(160) '%']));
end
%weather_absPressure = cell2mat(weather(9)) ;%surface pressure [mb]
weather_absPressure = weatherPressure(day_);

%weather_Temperature = cell2mat(weather(3));%surface temperature [deg C]
weather_Temperature = (weatherTemperature(day_) -32).*5/9;

%weather_VW = cell2mat(weather(5));%surface VW relative humidity [%]
weather_VW = weatherWV(day_);

%Delete duplicate time stamps in case weather station writes duplicates
[weather_dateTime,IA,~] = unique(weather_dateTime);
weather_Temperature = weather_Temperature(IA);
weather_absPressure = weather_absPressure(IA);
weather_VW = weather_VW(IA);

%interpolate surface data to analysis time grid
weather_Temperature = filloutliers(weather_Temperature,'linear');
weather_Temperature_interp = interp1(weather_dateTime,weather_Temperature,time_date_hms,'nearest','extrap');
weather_Temperature_interp = fillmissing(weather_Temperature_interp,'nearest');
weather_absPressure_interp = interp1(weather_dateTime,weather_absPressure,time_date_hms,'nearest','extrap');
weather_absPressure_interp = fillmissing(weather_absPressure_interp,'nearest');
weather_WV_interp = interp1(weather_dateTime,weather_VW,time_date_hms,'nearest','extrap');
weather_WV_interp = fillmissing(weather_WV_interp,'nearest');

%convert relative humidity to absolute humidity
% T = weather_Temperature_interp + 273.15;
% Pws = exp(77.3450+0.0057.*T-7235./T)./T.^8.2;
% weather_WV_interp = weather_WV_interp.*0.0022.*Pws./T./100;



end
