function [cohsig7] = RB_O2_770_PCA(T,P,nu_Range)
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

%Coherent s7
%load('RBPCA_1_15_21single.mat')
%spontaneous s6
load('RBPCA_11_2_21single.mat')

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

cohsig7I = ones(length(nu),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:))
        %cross_section(j,i,:) = interp1(nupermute,muY + M*thetapermute(:,i,j),nu_Rangepermute,'nearest'); 
        cohsig7I(:,i,j) = muY + M*thetapermute(:,i,j);  
        %cross_section(j,i,:) = muY + M*permute(theta(j,i,:),[3 2 1]); 
    end
end



cohsig7 = permute(cohsig7I(:,:,:),[3 2 1]);
normcohsig7 = trapz(cohsig7,3)*(nu_Range(2)-nu_Range(1))*100;
cohsig7 = cohsig7./normcohsig7;


end