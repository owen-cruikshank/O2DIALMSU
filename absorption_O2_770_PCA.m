function [absorption,cross_section] = absorption_O2_770_PCA(T,P,nu_Range,WV)
%File: absorption_O2_770_model_wavenumber.m
%Date: 02/28/2020
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -nu_Range:[1/cm] wavenumber scan. Dimentions: (1 x 1 x wavenumber)
%   -WV:[1/m^3] water vapor number density. Dimentions: (range x time)
%
%Outputs:
%   -absorption: [1/m] the atmopheric absorption of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time x wavenumber)
%   -sigma: [m^2] The voight absorption cross section of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time x wavenumber)
%   -f: [m] Absorption lineshape function, dimensions (range x time x wavnumber)


% % lambdaCenter = 769.7958;
% % lambdaWidth = 0.2;
% % lambdaDelta = .00002;
% % lambda = lambdaCenter-lambdaWidth/2:lambdaDelta:lambdaCenter+lambdaWidth/2;%[nm]
% % nu = 10^7./lambda;
% % 
% % 
% % %nu = permute(nu_scan,[3 1 2]);
% % nu = permute(nu,[3 1 2]);


%load variables needed
%load('PCA2.mat');
%load('PCA_5_18_20.mat');
%%%load('PCA_1_11_21.mat');
if length(nu_Range)==1
    load('PCA_1_11_21singleOnline.mat');
else
%load('PCA_1_11_21double.mat');
load('PCA_1_11_21single.mat');
end

% lambdaCenter = 769.7958;
% lambdaWidth = 0.2;
% lambdaDelta = .00002;
% lambda = lambdaCenter-lambdaWidth/2:lambdaDelta:lambdaCenter+lambdaWidth/2;%[nm]
% nu = 10^7./lambda;
% nu = permute(nu,[3 2 1]);

%cross_section = ones(length(T(:,1)),length(T(1,:)),length(nu_Range));
%cross_sectionI = ones(length(T(:,1)),length(T(1,:)),length(nu)));


%order
No = 20;
        normT = (T-muT)/sigmaT;
        normP = (P-muP)/sigmaP;
        
inc=1;
for n = 1:No
    for m = 1:n
        theta(:,:,inc) = normT(:,:).^(n-m) .* normP(:,:).^(m-1);
        inc=inc+1;
    end
end

%nupermute = permute(nu,[3 2 1]);
%nu_Rangepermute = permute(nu_Range,[3 2 1]);
thetapermute = permute(theta(:,:,:),[3 2 1]);

cross_sectionI = ones(length(nu),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:))
        %cross_section(j,i,:) = interp1(nupermute,muY + M*thetapermute(:,i,j),nu_Rangepermute,'nearest'); 
        cross_sectionI(:,i,j) = muY + M*thetapermute(:,i,j);  
        %cross_section(j,i,:) = muY + M*permute(theta(j,i,:),[3 2 1]); 
    end
end


if length(nu)==1
    cross_section = permute(cross_sectionI(:,:,:),[3 2 1]);
else
    %cross_section = permute(cross_sectionI(1:2:end,:,:),[3 2 1]);
cross_section = permute(cross_sectionI(:,:,:),[3 2 1]);
end
%rm = 1:length(T(:,1));
%ts = 1:length(T(1,:));

%cross_section = interp3(rm,ts,nupermute,cross_sectionI,rm,ts,nu_Rangepermute);

kB = 1.38065E-23;                       %[J/K] Boltzman's constant 
q_O2 = .2095;                           %[unitless] O2 atmospheric mixing ratio 
N_o2 = ((P*101325)./(kB*T)-WV) * q_O2;     %[molecule/m^3](t x r)O2 number density from atmopsheric number density and O2 mixing ratio
%N_o2 = ((P*101325)./(kB*T)) * q_O2;     %[molecule/m^3](t x r)O2 number density from atmopsheric number density and O2 mixing ratio

absorption = cross_section .* N_o2;     %[1/m](t x r x nu)absorption coefficeint of oxygen in the atmosphere
end