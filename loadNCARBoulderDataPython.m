function [Data] = loadNCARBoulderDataPython(spanDays,path)
yearB = num2str(year(spanDays));

fullDataPath = strcat(path,'mpd05.',yearB(3:4),sprintf('%02d',month(spanDays)),sprintf('%02d',day(spanDays)),'.Python.nc');

Data.Backscatter_Ratio = ncread(fullDataPath,'Backscatter_Ratio');%BSR unitless
Data.Backscatter_Ratio_mask = ncread(fullDataPath,'Backscatter_Ratio_mask');%True-masked, false-notmasked
Data.Aerosol_Backscatter_Coefficient = ncread(fullDataPath,'Aerosol_Backscatter_Coefficient');%'m^(-1) sr^(-1)'
%Data.Molecualr_Backscatter_Coefficient = Data.Aerosol_Backscatter_Coefficient/(Data.Backscatter_Ratio-1);%'m^(-1) sr^(-1)'
Data.HSRLMolecular_RayleighBrillioun = ncread(fullDataPath,'HSRLMolecular_RayleighBrillioun');%integral (sum=1) normalized Rayleigh-Brillioun spectrum.
Data.r_freqency = ncread(fullDataPath,'r_frequency');%Relative (offset) frequency axis in Hz

Data.HSRLCombined_scan_wavelength = ncread(fullDataPath,'HSRLCombined_scan_wavelength');
Data.HSRLCombined_TransmissionNorm = ncread(fullDataPath,'HSRLCombined_TransmissionNorm');
Data.HSRLMolecular_scan_wavelength = ncread(fullDataPath,'HSRLMolecular_scan_wavelength');
Data.HSRLMolecular_TransmissionNorm = ncread(fullDataPath,'HSRLMolecular_TransmissionNorm');

Data.HSRLMolecular_RangeFilterWidth = ncread(fullDataPath,'HSRLMolecular_RangeFilterWidth');
Data.HSRLMolecular_TimeFilterWidth = ncread(fullDataPath,'HSRLMolecular_TimeFilterWidth');

Data.Absolute_Humidity = ncread(fullDataPath,'Absolute_Humidity');
Data.Absolute_Humidity_mask = ncread(fullDataPath,'Absolute_Humidity_mask');


Data.Temperature_Model = ncread(fullDataPath,'Temperature_Model');
Data.Pressure_Model = ncread(fullDataPath,'Pressure_Model');

Data.Surface_AbsHum = ncread(fullDataPath,'Surface_AbsHum');
Data.Surface_Pressure = ncread(fullDataPath,'Surface_Pressure');
Data.Surface_Temperature = ncread(fullDataPath,'Surface_Temperature');

Data.range = ncread(fullDataPath,'range');
Data.time = ncread(fullDataPath,'time');

Data.HSRLMolecular_Wavelength_median = ncread(fullDataPath,'HSRLMolecular_Wavelength_median');
Data.HSRLCombined_Wavelength_median = ncread(fullDataPath,'HSRLCombined_Wavelength_median');



end