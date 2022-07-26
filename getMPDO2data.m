% This function will extract the desired MPD #5 data (useful if multiple days
% are needed)
function [o2on,o2onr,o2ont,o2off,o2offr,o2offt] = getMPDO2data(filename)

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

 o2on = O2_Online_Backscatter_Channel_Raw_Data;
 o2onr = range_O2_Online_Backscatter_Channel_Raw_Data;
 o2ont = time_O2_Online_Backscatter_Channel_Raw_Data;
 
 o2off = O2_Offline_Backscatter_Channel_Raw_Data;
 o2offr = range_O2_Offline_Backscatter_Channel_Raw_Data;
 o2offt = time_O2_Offline_Backscatter_Channel_Raw_Data;
