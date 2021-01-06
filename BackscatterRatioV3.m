function [Bm,Ba,BR]=...
    BackscatterRatioV3(time,range,Nc,Nm,T,P,Wavelength)
%Inputs
%time is constant time grid that all data is interpolated or integrated to
%range in km
%Nc and Nm are counts for combined and molecular channels
%T is temperature in Kelvin
%P is pressure in atmospheres
%Wavelength is wavelength vector from laser locking, O2Offline

%Outputs
%Bm is molecular backscatter coefficient
%Ba is aerosol backscatter coefficient
%BR is backscatter ratio




%% Molecular Backscatter Coefficient
Bm=zeros(size(Nc));
Bm=374280*P./(T.*Wavelength.^4);
Bm(end,:)=Bm(end-1,:);


%% Cmm
%I found that Cmm is highly temperature dependent.
%I used a polynomial fit since interpolating took too long
%before 8/27/2020
% CmmPolyfit=load('Cmm(T).mat'); 
% CmmPolynomial=load('CmmPolynomial.mat');
%after 8/27/2020
load('CalibrationData.mat')
Cmm=polyval(CmmPolynomialConstants,T);


%Calibration Constants
Cac=1; %Aerosol in Combined
Cmc=0.998; %Molecular in Combined
Cmc=0.977;
Cmc=0.965;
Cmc=0.975;
Cam=0.00056; %Aerosol in Molecular
Cam=0.00095;
Cam=2.2e-3;
Cam=3.94e-4;

%Beamsplitter and Detector Efficiencies
eta_c=.5217;
eta_c=.4506;
eta_c=.397;
eta_c=.406;
eta_m=.4783;
eta_m=.5494;
eta_m=0.603;
eta_m=0.594;
%Aerosol Backscatter Coefficient
Ba=Bm.*(Cmm*eta_m.*(Nc)./(Nm)-Cmc*eta_c)./(Cac*eta_c-Cam*eta_m*(Nc)./(Nm));

%Backscatter Ratio
BR=(Ba./Bm)+1;


%Smoothing the Backscatter Ratio
% BR=smoothdata(BR,1,'g',20);

end