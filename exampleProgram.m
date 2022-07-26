%example program
clear

%===== Varibles
%=Date
date_start = datetime(2021,3,28,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2021,4,3,'TimeZone','UTC');%yyyy,mm,dd

%date_start = datetime(2021,4,28,'TimeZone','UTC');%yyyy,mm,dd
%date_end = datetime(2021,5,3,'TimeZone','UTC');%yyyy,mm,dd

%date_start = datetime(2021,5,3,'TimeZone','UTC');%yyyy,mm,dd
%date_end = datetime(2021,5,5,'TimeZone','UTC');%yyyy,mm,dd

span_days = date_start:date_end;

%=Time and range
Options.intTime = 10;  %[min] Integration time
Options.intRange = 1; %[bins] Integration range

Options.t_avg = 1;     %[bins] Time smoothing bins
Options.oversample = 1; %[bins] Range smoothing bins



%====================
%==== Constants =====
%====================
Constant.g0 = 9.80665;       %[m/s^2] Gravitational acceleration 
Constant.M_air = 0.0289644;  %[kg/mol] Molar mass of Earth's air 
Constant.R = 8.3144598;      %[J/(mol*K)] Universal gas constant 

Constant.c = 2.99792458E8;           %[m/s] Speed of light 
Constant.kb = 1.38065E-23;           %[J/K][m^2 kg s-2 K-1] Boltzman's constant 
Constant.h = 6.626E-34;              %[Js] Planck's constant 
Constant.mo2 = 5.314E-26;            %[kg] Mass O2 molecule 
Constant.mWV = 2.9915e-26;           %[kg] Mass H2O molecule
Constant.m_air = 4.792E-26;          %[kg] Mass of air
Constant.q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
Constant.No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)

[Range,Time,Counts,Sonde,Model,Spectrum,BSR,Data] = loadMSUdata(span_days,Options,Constant);