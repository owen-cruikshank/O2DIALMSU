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
%     CmmPolynomial=load('CmmPolynomial.mat'); 
%     LidarData.Cmm=polyval(CmmPolynomial.Constants,WeatherData.Temperature);
    file=pwd;
    cd("F:\Research\Calibration_Data")
    load('Rayleigh1229.mat');
    cd(file)
    
    %Constants
    LidarData.Cam=Results.Cam(1);
    LidarData.Cac=Results.Cac(1);

    %Polynomial Fit
     LidarData.Cmm=(Results.PolynomialfitCmm(1)).*...
    (WeatherData.Temperature).^4+(Results.PolynomialfitCmm(2)).*...
    (WeatherData.Temperature).^3+(Results.PolynomialfitCmm(3)).*...
    (WeatherData.Temperature).^2+(Results.PolynomialfitCmm(4)).*...
    (WeatherData.Temperature)+(Results.PolynomialfitCmm(5)); 
    %Aerosol in Combined
     LidarData.Cmc=(Results.PolynomialfitCmc(1)).*...
    (WeatherData.Temperature).^4+(Results.PolynomialfitCmc(2)).*...
    (WeatherData.Temperature).^3+(Results.PolynomialfitCmc(3)).*...
    (WeatherData.Temperature).^2+(Results.PolynomialfitCmc(4)).*...
    (WeatherData.Temperature)+(Results.PolynomialfitCmc(5));    
    %Aerosol in Molecular

%     %Beamsplitter and Detector Efficiencies
     LidarData.EtaCombined=Results.CombinedEfficiency;
     LidarData.EtaMolecular=Results.MolecularEfficiency;

    %Aerosol Backscatter Coefficient
    LidarData.AerosolBackscatterCoefficient=LidarData.MolecularBackscatterCoefficient.*(LidarData.Cmm*LidarData.EtaMolecular.*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts)-LidarData.Cmc*LidarData.EtaCombined)./(LidarData.Cac*LidarData.EtaCombined-LidarData.Cam*LidarData.EtaMolecular*(LidarData.OfflineCombinedAverageCounts)./(LidarData.OfflineMolecularAverageCounts));

    %Backscatter Ratio
    LidarData.BackscatterRatio=(LidarData.AerosolBackscatterCoefficient./LidarData.MolecularBackscatterCoefficient)+1;

end