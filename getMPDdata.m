% This function will extract the desired MPD data (useful if multiple days
% are needed)
%function [ba,bar,bat,T,Tr,Tt,P,Pr,Pt,WV,WVr,WVt] = getMPDdata(filename)       % No masks
function [ba,bam,bar,bat,T_surf,T,Tr,Tt,P_surf,P,Pr,Pt,WV_surf,WVm_surf,WV,WVm,WVr,WVt] = getMPDdata(filename) % Includes masks
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

 ba = Aerosol_Backscatter_Coefficient;
 bam = Aerosol_Backscatter_Coefficient_mask;
 bar = range_Aerosol_Backscatter_Coefficient;
 bat = time_Aerosol_Backscatter_Coefficient;
 
 T_surf = Surface_Temperature_HSRL;
 T = Temperature;
 Tr = range_Temperature;
 Tt = time_Temperature;
 
 P_surf = Surface_Pressure_HSRL;
 P = Pressure;
 Pr = range_Pressure;
 Pt = time_Pressure;
 
 WV_surf = Surface_Absolute_Humidity;
 WVm_surf = Surface_Absolute_Humidity_mask;
 WV = Absolute_Humidity;
 WVm = Absolute_Humidity_mask;
 WVr = range_Absolute_Humidity;
 WVt = time_Absolute_Humidity;