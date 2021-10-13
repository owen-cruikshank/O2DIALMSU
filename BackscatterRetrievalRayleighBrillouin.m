function [LidarData]=BackscatterRetrievalRayleighBrillouin(LidarData,WeatherData)
%inputs-
%Lidar Data with Corrected Counts
%Weather Data with Temperature and Pressure
%outputs-
%Lidar Data with Backscatter Coefficients and Backscatter Ratio

    %Create Molecular Backscatter Coefficient
    LidarData.MolecularBackscatterCoefficient=zeros(size(LidarData.OfflineCombinedAverageCounts));
    LidarData.MolecularBackscatterCoefficient=374280*WeatherData.Pressure./(WeatherData.Temperature.*770.1085.^4);
    LidarData.MolecularBackscatterCoefficient(end,:)=LidarData.MolecularBackscatterCoefficient(end-1,:);


    % Finding Cmm- The Constant of Molecular Scattering in the Molecular Channel
    %I found that Cmm is highly temperature dependent.
    %I used a polynomial fit
%     CmmPolynomial=load('CmmPolynomial.mat'); 
%     LidarData.Cmm=polyval(CmmPolynomial.Constants,WeatherData.Temperature);
    %file=pwd;
    %cd("F:\Research\Calibration_Data")
    load('CalibrationTables0309.mat');
    %cd(file)
    P=Results.Pressure*0.009869233; 
    T=Results.Temperature;
    Eta_m=Results.MolecularEfficiency;
    Eta_c=Results.CombinedEfficiency;
    Cmm=Results.Cmm;
    Cmc=Results.Cmc; 
    Cam=Results.Cam(1); 
    Cac=Results.Cac(1);
    
    
LidarData.Cmm=zeros(length(LidarData.Range),length(LidarData.Time));
LidarData.Cmc=zeros(length(LidarData.Range),length(LidarData.Time));

    for i=1:length(LidarData.Time)
       for j=1:length(LidarData.Range)
        LidarData.Cmm(j,i)=interp2(P,T,Cmm,WeatherData.Pressure(j,i),WeatherData.Temperature(j,i));
        LidarData.Cmc(j,i)=interp2(P,T,Cmc,WeatherData.Pressure(j,i),WeatherData.Temperature(j,i));
       end
    end

    %Calibration Constants
    LidarData.Cac=Cac; %Aerosol in Combined
    LidarData.Cam=Cam; %Aerosol in Molecular

    %Beamsplitter and Detector Efficiencies
    LidarData.EtaCombined=Eta_c;
    LidarData.EtaMolecular=Eta_m;

    %Aerosol Backscatter Coefficient
    LidarData.AerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*(LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts)-LidarData.Cmc*LidarData.EtaCombined) ...
    ./(LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts));

    LidarData.o2off_mol_corr = ((LidarData.OfflineMolecularAverageCounts ./Eta_m)-LidarData.Cam.*LidarData.OfflineCombinedAverageCounts)./(LidarData.Cmm-LidarData.Cam);
    %Backscatter Ratio
    LidarData.BackscatterRatio=(LidarData.AerosolBackscatterCoefficient./LidarData.MolecularBackscatterCoefficient)+1;
    %set values less than 1 equal to 1
    LidarData.BackscatterRatio(LidarData.BackscatterRatio<1)=1;
end