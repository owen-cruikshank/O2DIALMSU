% This function will extract the desired MPD data (useful if multiple days
% are needed)
function [T, P, r_TP,t_single_SGP] = getSGPdata(filename)
%function [ba,bam,bar,bat,T,Tr,Tt,P,Pr,Pt,WV,WVm,WVr,WVt] = getMPDdata(filename)
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

T = tdry + 273.15;  % Radiosonde temperature Celsius to [Kelvin]
P = pres./1013.25;  % Radiosonde pressure mbar to [atm]
r_TP = alt;         % Range [m]
% Datetime of a single SGP radiosonde
t_single_SGP = datetime((double(base_time) + time),'convertfrom','posixtime');
