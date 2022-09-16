function absorption = absorption_O2_770_model_wavenumber(T,P,nu_Range,WV)
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

c = 2.99792458E8;                       %[m/s] speed of light 
kB = 1.38065E-23;                       %[J/K] Boltzman's constant 
h = 6.62607004E-34;                     %[Js] Planck's constant

q_O2 = .2095;                           %[unitless] O2 atmospheric mixing ratio 

mo2 = 5.313525282632483E-26;            % Mass O2 molecule [kg]

N_o2 = ((P*101325)./(kB*T)-WV) * q_O2;     %[molecule/m^3](t x r)O2 number density from atmopsheric number density and O2 mixing ratio
%N_o2 = ((P*101325)./(kB*T)) * q_O2;     %[molecule/m^3](t x r)O2 number density from atmopsheric number density and O2 mixing ratio

%reference T and P
T0 = 296;                               %[K]
P0 = 1.0;                               %[atm]

parameters = fopen(fullfile('CalibrationData','O2_line_parameters.out'),'r');   %open file containing HITRAN information
fmt = ['%1d %1d %f %e %e %f %f %f %f %f'];          %format of HITRAN file
O2_parameters =fscanf(parameters,fmt,[10 inf]);     %place HITRAN parameters in vector a
fclose(parameters);                                 %close file
O2_parameters = O2_parameters';                     %transpose matrix to correct format

[rL, tL] = size(T);                                 %length of range vector x length of time vector

%lineshape = zeros(rL,tL);                           %preallocate matricies
%Voight_profile = zeros(rL,tL);
%cross_section = zeros(rL,tL);
absorption = zeros(rL,tL,length(nu_Range));%preallocate 
integralV = zeros(rL,tL,length(nu_Range));%preallocate 

nu_Range = nu_Range * 100;                          %change nu_Range from [1/cm] to [1/m]

strength_threshold = 1*10^(-27);                    %[cm / molecule] line strength threshold
strength_threshold = 0;                    %[cm / molecule] line strength threshold

t = -10:.2:10;                                      %Relative freqency to integrate over
t_4D = permute(t,[4 3 1 2]);                        %shift t to put it in third dimestion
figure
hold on

for i = 1:size(O2_parameters,1)                     %loop over all line parameters

    nu_O2 = O2_parameters(i,3);                     %[1/cm]
    nu_O2 = nu_O2 * 100;                            %[1/m] absoption wavenumber
    S0_O2 = O2_parameters(i,4);                     %[cm/molecule] line strength
    S0_O2 = S0_O2 / 100;                            %[m/molecule] absoption line strength at T=296K
    %S0_O2 = S0_O2 * 0.75;
    gamma_L = O2_parameters(i,6);                   %[1/cm] linewidth
    gamma_L = gamma_L * 100;                        %[1/m] Air (Lorentz) broadened linewidth
    n_air = O2_parameters(i,9);                     %[unitless] linewidth temperature dependence
    E_lower = O2_parameters(i,8);                   %[1/cm] ground state energy
    E_lower = E_lower * 100;                        %[1/m] ground state energy
    gamma_D = O2_parameters(i,7);                   %[1/cm]self (Doppler) broadened linewidth
    gamma_D = gamma_D * 100;                        %[1/m]
    delta_air = O2_parameters(i,10);                %[1/cm/atm] pressure shift induce by air, at p=1atm
    %delta_air = -0.0093; %override
    %delta_air = -0.013; %override
    %delta_air = -0.02; %override
    delta_air = delta_air * 100;                    %[1/m/atm]
    
    %if S0_O2 * 100 > strength_threshold             %Do not compute cross section for line if it id below the threshold
    %if nu_O2 >= nu_Range(:,:,1) && nu_O2 <= nu_Range(:,:,end) && S0_O2 * 100 > strength_threshold %Do not compute if line is out of specified wavnumber range

        nuShifted = nu_O2 + delta_air .* P;         %[1/m] (r x t) Shift line center based on atmopheric pressure
        %nuShifted = nu_O2 ;%+ delta_air .* P;         %[1/m] (r x t) Shift line center based on atmopheric pressure

        %temperature shifted line strength
        a_o2 = S0_O2.*(T0./T);
        c_o2 = exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower);
        ST_O2 = a_o2.*c_o2;                                     %[m/molecule](t x r) O2 line strength adjusted for temperature shift from T0
        gamma_L_T = gamma_L * (P/P0).*((T0./T).^n_air);         %[1/m](t x r) Lorentz lineshape adjusted for temperature and pressure shift
        gamma_D_T = (nuShifted/c).*sqrt(2*kB*T*log(2)/mo2); %[1/m](t x r) Dopper lineshape due to temperature
        
        %voight lineshape
        x = ((nu_Range-nuShifted)./gamma_D_T) * sqrt(log(2));   %[none](t x r x nu)
        y = (gamma_L_T./gamma_D_T) * sqrt(log(2));              %[none](t x r)
        K = (ST_O2./gamma_D_T) * sqrt(log(2)/pi);               %[m^2 / molecule (t x r)



        %M = 4;                                                 % Number of parallel workers
        %parfor (time_i = 1:tL, M)                              % Loop over time vector
        parfor time_i = 1:tL                                       % Loop over time vector
            %integration_function = exp(-t_4D.^2)./(y(:,time_i).^2 + (x(:,time_i,:)-t_4D).^2);   %[none] create integration function to integrate over
            integralV(:,time_i,:) = trapz(t,exp(-t_4D.^2)./(y(:,time_i).^2 + (x(:,time_i,:)-t_4D).^2),4);                            %[none] integrate over t
        end

        %f = log(2).*pi^(-3/2).*gamma_L_T./gamma_D_T.^2.*integralV;                  %[m](t x r x nu) Absorption lineshape 

        %sigma =   ;                                                      %[m^2](t x r x nu) Voight absorption cross section


        %lineshape = f + lineshape;                          %add on to previous lineshape
        %cross_section = sigma + cross_section;              %add on to previous cross_section
        line = (y/pi).*integralV .* K .* N_o2;
        absorption = absorption + (y/pi).*integralV .* K .* N_o2; %[1/m](t x r x nu)absorption coefficeint of oxygen in the atmosphere
        semilogy(permute(nu_Range,[3 2 1]),permute(line(1,1,:),[3 2 1]))
   % end
    
end 
hold off
end