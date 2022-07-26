%weather example
span_days = datetime(2020,2,22,'TimeZone','UTC');
t_end = 86340; %[s] Ending time  
t_step = 60;                               %[s] Time step
ts = 0:t_step:t_end;                       %[s] Time vector
path = 'D:\Owen\OneDrive - Montana State University - Bozeman\research s19\o2DIAL_data';
[weather_Temperature_interp, weather_absPressure_interp, weather_VW_interp] = ORSLweather(span_days,ts,path)