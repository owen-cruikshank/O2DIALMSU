% This function will extract the desired SGP tower data
function [T_02m, T_25m, T_60m, P_02m, P_25m, P_60m, t_SGPTWR] = getSGPTWRdata(filename)

% ========================================================================
% ======== Import NetCDF data ========
% ========================================================================

f = ncinfo(filename);
nvars = length(f.Variables);
for k = 1:nvars
   varname=f.Variables(k).Name;
   %disp(['Reading:  ' varname]);
   eval([varname ' = ncread(''' filename ''',''' varname ''');']);
   % Suppress warning that names are limited to 63 characters
   warning('off', 'MATLAB:namelengthmaxexceeded')
end

T_02m = temp_02m + 273.15;  % Tower temperature at 2 m (Celsius to [Kelvin])
T_25m = temp_25m + 273.15;  % Tower temperature at 25 m (Celsius to [Kelvin])
T_60m = temp_60m + 273.15;  % Tower temperature at 60 m (Celsius to [Kelvin])
P_02m = pres_02m./1013.25;  % Tower pressure at 2 m (mbar to [atm])
P_25m = pres_25m./1013.25;  % Tower pressure at 25 m (mbar to [atm])
P_60m = pres_60m./1013.25;  % Tower pressure at 60 m (mbar to [atm])
% Datetime of a single SGP radiosonde
t_SGPTWR = datetime((double(base_time) + time),'convertfrom','posixtime');
