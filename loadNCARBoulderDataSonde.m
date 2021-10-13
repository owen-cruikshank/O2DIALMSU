function [Data] = loadNCARBoulderDataSonde(spanDays,path)
yearB = num2str(year(spanDays));

fullDataPath = strcat(path,'Marshall_Field_Site_',yearB,sprintf('%02d',month(spanDays)),sprintf('%02d',day(spanDays)),'_161051.nc');
topInd = 2800;
try
    %refAlt = ncread(fullDataPath,'reference_alt');
    range = ncread(fullDataPath,'alt');
    range = range(1:topInd);
    Data.rm = {range-range(1)};
    temperature = ncread(fullDataPath,'tdry')+273.15;%C
    temperature = temperature(1:topInd);
    Data.temperature = {temperature};
    pressure = ncread(fullDataPath,'pres')/1013.2501;
    pressure = pressure(1:topInd);
    Data.pressure = {pressure};%hPa
    ts = ncread(fullDataPath,'time');%seconds since 2020-09-02T16:10:51 UTC
    ts = ts(1:topInd);
    %Data.rh = ncread(fullDataPath,'rh');%relative humidity %
    %date = datetime(2020,09,02,16,10,51,'TimeZone','UTC');
    date = datetime(2020,09,02,16,10,51);
    Data.date = {date};
    Data.date_ts = {date+seconds(ts)};
catch
    Data.temperature = {[]};
    Data.rm = {[]};
    Data.pressure = {[]};
    Data.date_ts = {[]};
end

end