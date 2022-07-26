function [Bm,Ba,BR,T,P]=...
    BackscatterRatio(time,dt,range,dr,Nc,Nm,GroundPressure,GroundTemperature,Wavelength)
%Inputs
%time is constant time grid that everythin is interpolated or integrated to
%dt is time size, I used minutes
%range and dr in km
%GroundPressure and GroundTemperature from weather station
%Wavelength is wavelength vector from laser locking, O2Offline
%Outputs
%Bm is molecular backscatter coefficient
%Ba is aerosol backscatter coefficient
%BR is backscatter ratio
%T is temperature estimation used
%P is pressure estimation used
%

%% Temperature and Pressure models
T=zeros(size(Nc));
P=zeros(size(Nc));

%Temperature estimated to be linear until 12 km then constant
lapserate=9.8;
lapserate=6.5;

T(:,1)=GroundTemperature+273.15;
troposphereindex=find(range<=12,1,'last');
T(:,2:troposphereindex)=T(:,1)-lapserate*range(2:troposphereindex);

for i=troposphereindex:length(range)-1
    T(:,i)=T(:,troposphereindex);
end

%Pressure from exponential atmospheric pressure model
P(:,1)=GroundPressure(:)/10;
P(:,2:end)=P(:,1).*(T(:,1)./T(:,2:end)).^-5.2199;

%% Molecular Backscatter Coefficient
Bm=zeros(size(Nc));
for i=1:length(time)
    Bm(i,:)=374280*P(i,:)./T(i,:)*(Wavelength(i))^(-4);
end

Bm(:,end)=Bm(:,end-1);

%% Aerosol Backscatter Coefficients

%Callibration Coefficients
Cac=1; %Aerosol in Combined
Cmc=0.92; %Molecular in Combined
Cmm=0.2; %Molecular in Molecular
Cam=0.0005; %Aerosol in Molecular

Cac=1; %Aerosol in Combined
Cmc=0.995; %Molecular in Combined
Cmm=0.18; %Molecular in Molecular
Cam=0.0017; %Aerosol in Molecular

%Aerosol Backscatter Coefficient
Ba=Bm.*(Cmm*(Nc)./(Nm)-Cmc)./(Cac-Cam*(Nc)./(Nm));

BR=(Ba./Bm)+1;

end