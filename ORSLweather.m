function [weather_Temperature_interp, weather_absPressure_interp, weather_VW_interp] = ORSLweather(span_days,ts,path)
%File: ORSLweather.m
%Date: 04/19/2020
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

%Seperate datetime vector to month and year
month_span = month(span_days);
if month_span < 10
    month_span = cat(2,'0',num2str(month_span));
else
    month_span = num2str(month_span);
end
year_span = year(span_days);
year_span = num2str(year_span);

%Read data in from file using month and year
file_date = cat(2,month_span,year_span,'.csv');
weatherFilePath = fullfile(path,'Weather station data',file_date);
weatherFile = fopen(weatherFilePath,'r');
formatSpec='%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
weather = textscan(weatherFile,formatSpec,'Delimiter',',','HeaderLines',1);
fclose(weatherFile);
%create datetime vector for date from weather
weather_dateTime =  strcat(weather{1},weather{2});
weather_dateTime = datetime(weather_dateTime,'TimeZone','UTC','InputFormat','MM/dd/yyyyHH:mm','Format','d-MMM-y HH:mm:ss Z');
%convert Denver to UTC
%weather_dateTime.TimeZone = 'UTC';

%find data we want and convert to matricies
day_ = day(weather_dateTime) == day(span_days);
day_ = logical(sum(day_,2));
weather_dateTime = weather_dateTime(day_);
weather_absPressure = cell2mat(weather(5));%surface pressure [mbar]
weather_absPressure = weather_absPressure(day_);
weather_Temperature = cell2mat(weather(7));%surface temperature [deg C]
weather_Temperature = weather_Temperature(day_);
weather_VW = cell2mat(weather(13));%surface VW number density [# cm-3]
weather_VW = weather_VW(day_);

%Delete duplicate time stamps in case weather station writes duplicates
[weather_dateTime,IA,IC] = unique(weather_dateTime);
weather_Temperature = weather_Temperature(IA);
weather_absPressure = weather_absPressure(IA);
weather_VW = weather_VW(IA);

%Filloutliers
%weather_Temperature = filloutliers(weather_Temperature,'linear');
weather_Temperature = filloutliers(weather_Temperature,'nearest','movmedian',50);

%interpolate surface data to analysis time grid
weather_Temperature_interp = interp1(weather_dateTime,weather_Temperature,time_date_hms,'nearest','extrap');
weather_absPressure_interp = interp1(weather_dateTime,weather_absPressure,time_date_hms,'nearest','extrap');
weather_VW_interp = interp1(weather_dateTime,weather_VW,time_date_hms,'nearest','extrap');

% weather_Temperature_interp = interp1(weather_dateTime,weather_Temperature,time_date_hms,'nearest');
% weather_absPressure_interp = interp1(weather_dateTime,weather_absPressure,time_date_hms,'nearest');
% weather_VW_interp = interp1(weather_dateTime,weather_VW,time_date_hms,'nearest');
weather_VW_interp=0;
end