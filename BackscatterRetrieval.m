function [LidarData]=BackscatterRetrieval(LidarData,WeatherData)
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
    CmmPolynomial=load('CmmPolynomial.mat'); 
    LidarData.Cmm=polyval(CmmPolynomial.Constants,WeatherData.Temperature);


    %Calibration Constants
    LidarData.Cac=1; %Aerosol in Combined
    LidarData.Cmc=0.977; %Molecular in Combined
    LidarData.Cam=0.00095; %Aerosol in Molecular

    %Beamsplitter and Detector Efficiencies
    LidarData.EtaCombined=.397;
    LidarData.EtaMolecular=.603;

    %Aerosol Backscatter Coefficient
    LidarData.AerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*(LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts)-LidarData.Cmc*LidarData.EtaCombined)./(LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts));

    %Backscatter Ratio
    LidarData.BackscatterRatio=(LidarData.AerosolBackscatterCoefficient./LidarData.MolecularBackscatterCoefficient)+1;

end