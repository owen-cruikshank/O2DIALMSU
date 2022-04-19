function [LidarData]=BackscatterRetrievalBoulder062021(LidarData,WeatherData)
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
    %%%%%load('CalibrationTablesBoulder062021.mat');
    load('CalibrationTablesBoulderSponS6062021.mat');
    BoulderHSRLcoefficents062021 = BoulderHSRLcoefficentsSponS6_062021;
    %cd(file)
    P=BoulderHSRLcoefficents062021.P; 
    T=permute(BoulderHSRLcoefficents062021.T,[2 1]);
    Eta_m=1;
    Eta_c=1;
    Cmm=permute(BoulderHSRLcoefficents062021.Cmm,[3 2 1]);
    Cmc=permute(BoulderHSRLcoefficents062021.Cmc,[3 2 1]); 
    Cam=BoulderHSRLcoefficents062021.Cam; 
    Cac=BoulderHSRLcoefficents062021.Cac;
    
    
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

    %wavelength corrected for wv
    LidarData.MolecularBackscatterCoefficient828 = LidarData.MolecularBackscatterCoefficient*770^4/828^4;
    
    %Aerosol Backscatter Coefficient
    LidarData.UnmaskedAerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*...
        (LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.OfflineCombinedTotalCounts)./...
        (LidarData.OfflineMolecularTotalCounts)-LidarData.Cmc*LidarData.EtaCombined)./...
        (LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*...
        (LidarData.OfflineCombinedTotalCounts)./(LidarData.OfflineMolecularTotalCounts));

    LidarData.UnmaskedAerosolBackscatterCoefficient828 = LidarData.UnmaskedAerosolBackscatterCoefficient*770/828;

    %Backscatter Ratio
    LidarData.UnmaskedBackscatterRatio=(LidarData.UnmaskedAerosolBackscatterCoefficient./...
        LidarData.MolecularBackscatterCoefficient)+1;

        LidarData.UnmaskedBackscatterRatio828=(LidarData.UnmaskedAerosolBackscatterCoefficient828./...
        LidarData.MolecularBackscatterCoefficient828)+1;

    
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