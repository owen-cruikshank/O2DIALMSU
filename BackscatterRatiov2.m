function [Bm,Ba,BR]=...
    BackscatterRatiov2(time,range,Nc,Nm,T,P,Wavelength)
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
% Bm=zeros(size(Nc));
% for i=1:length(time)
%     Bm(:,i)=374280*P(:,i)./(T(:,1)*(Wavelength(i))^(4));
% end

Bm=374280*P./(T.*Wavelength.^4);


Bm(end,:)=Bm(end-1,:);


%% Cmm
%I found that Cmm is highly temperature dependent.
%I used a polynomial fit since interpolating took too long
CmmPolyfit=load('CmmPolyfit.mat'); 
Cmm=polyval(CmmPolyfit.CmmPolyfit,T);


%Calibration Constants
Cac=1; %Aerosol in Combined
Cmc=0.995; %Molecular in Combined
Cam=0.0017; %Aerosol in Molecular

%Beamsplitter Efficiencies
eta_c=.3896;
eta_m=.6104;

%Aerosol Backscatter Coefficient
Ba=Bm.*(Cmm*eta_m.*(Nc)./(Nm)-Cmc*eta_c)./(Cac*eta_c-Cam*eta_m*(Nc)./(Nm));

%Backscatter Ratio
BR=(Ba./Bm)+1;


%Smoothing the Backscatter Ratio
%BR=smoothdata(BR,1,'g',10);

end