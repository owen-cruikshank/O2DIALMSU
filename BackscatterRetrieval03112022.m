function [LidarData]=BackscatterRetrieval03112022(LidarData,WeatherData)
%inputs-
%Lidar Data with Corrected Counts
%Weather Data with Temperature and Pressure
%outputs-
%Lidar Data with Backscatter Coefficients and Backscatter Ratio
%% Load in calibration data
    % Finding Cmm- The Constant of Molecular Scattering in the Molecular Channel
    %I found that Cmm is highly temperature dependent.
    %I used a polynomial fit
%     CmmPolynomial=load('CmmPolynomial.mat'); 
%     LidarData.Cmm=polyval(CmmPolynomial.Constants,WeatherData.Temperature);
    %file=pwd;
    %cd("F:\Research\Calibration_Data")
   % load(fullfile('CalibrationData','CalibrationTables0702.mat'));
    load(fullfile('CalibrationData','CalibrationScan03112022.mat'));
    %cd(file)
    P=Results.Pressure*0.009869233; 
    T=Results.Temperature;
    Eta_m=Results.eta_m;
    Eta_c=Results.eta_c;
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
    
%% Unmasked Version
    %Create Molecular Backscatter Coefficient
    LidarData.MolecularBackscatterCoefficient=9.94266e-7*(WeatherData.Pressure)./(WeatherData.Temperature);

    
    %Aerosol Backscatter Coefficient
    LidarData.UnmaskedAerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*...
        (LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.OfflineCombinedTotalCounts)./...
        (LidarData.OfflineMolecularTotalCounts)-LidarData.Cmc*LidarData.EtaCombined)./...
        (LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*...
        (LidarData.OfflineCombinedTotalCounts)./(LidarData.OfflineMolecularTotalCounts));

    %Backscatter Ratio
    LidarData.UnmaskedBackscatterRatio=(LidarData.UnmaskedAerosolBackscatterCoefficient./...
        LidarData.MolecularBackscatterCoefficient)+1;

    
%% Masked Version
%     %Aerosol Backscatter Coefficient
%     LidarData.MaskedAerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*...
%         (LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.MaskedOfflineCombinedTotalCounts)./...
%         (LidarData.MaskedOfflineMolecularTotalCounts)-LidarData.Cmc*LidarData.EtaCombined)./...
%         (LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*...
%         (LidarData.MaskedOfflineCombinedTotalCounts)./(LidarData.MaskedOfflineMolecularTotalCounts));
% 
%     %Backscatter Ratio
%     LidarData.MaskedBackscatterRatio=(LidarData.MaskedAerosolBackscatterCoefficient./...
%         LidarData.MolecularBackscatterCoefficient)+1;

    
    
% %     fields={'Cmm','Cmc','Cac','Cam','EtaCombined','EtaMolecular',...
% %         'MaskedOnlineMolecularTotalCounts','MaskedOnlineCombinedTotalCounts',...
% %         'MaskedOfflineMolecularTotalCounts','MaskedOfflineCombinedTotalCounts',...
% %         'OnlineMolecularTotalCounts','OfflineMolecularTotalCounts',...
% %         'OfflineCombinedTotalCounts','OnlineCombinedTotalCounts'};
%         fields={'Cmm','Cmc','Cac','Cam','EtaCombined','EtaMolecular',...
%         'OnlineMolecularTotalCounts','OfflineMolecularTotalCounts'};
%         
%     LidarData=rmfield(LidarData,fields);
end