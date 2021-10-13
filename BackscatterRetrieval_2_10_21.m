function [LidarData]=BackscatterRetrieval_2_10_21(LidarData,WeatherData)
%inputs-
%Lidar Data with Corrected Counts
%Weather Data with Temperature and Pressure
%outputs-
%Lidar Data with Backscatter Coefficients and Backscatter Ratio

    %Create Molecular Backscatter Coefficient
    LidarData.MolecularBackscatterCoefficient=zeros(size(LidarData.OfflineCombinedAverageCounts));
    LidarData.MolecularBackscatterCoefficient=374280*WeatherData.Pressure./(WeatherData.Temperature.*770.1085.^4);
    LidarData.MolecularBackscatterCoefficient(end,:)=LidarData.MolecularBackscatterCoefficient(end-1,:);

    load('Calibration0203.mat')

    % Finding Cmm- The Constant of Molecular Scattering in the Molecular Channel
    %I found that Cmm is highly temperature dependent.
    %I used a polynomial fit
    %CmmPolynomial=load('CmmPolynomial.mat'); 
    %LidarData.Cmm=polyval(CmmPolynomial.Constants,WeatherData.Temperature);
    LidarData.Cmm=polyval(Results.PolynomialfitCmm,WeatherData.Temperature);


    %Calibration Constants
    LidarData.Cac=Results.Cac(1); %Aerosol in Combined
    LidarData.Cmc=Results.Cmc(1); %Molecular in Combined
    LidarData.Cam=Results.Cam(1); %Aerosol in Molecular

    %Beamsplitter and Detector Efficiencies
    LidarData.EtaCombined=Results.CombinedEfficiency;
    LidarData.EtaMolecular=Results.MolecularEfficiency;

    %Aerosol Backscatter Coefficient
    LidarData.AerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*(LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts)-LidarData.Cmc*LidarData.EtaCombined)./(LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts));

    %Backscatter Ratio
    LidarData.BackscatterRatio=(LidarData.AerosolBackscatterCoefficient./LidarData.MolecularBackscatterCoefficient)+1;

end