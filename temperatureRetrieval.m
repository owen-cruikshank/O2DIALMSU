function [T_final,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion] =  temperatureRetrieval(T,ts,rm,P,WV,nu_scan,alpha_O2,SNRm,cloud_SDm_above)
%File: temperatureRetrieval.m
%Date: 03/16/2020
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -ts:[s] (1 x time) vector of time dimension
%   -rm:[m] (range x 1) vector of range dimension
%   -WV:[molecules/m^3] scalar or (range x time) vector of water vapor
%   number density
%   -nu_scan:[1/cm] wavenumber of online laser. Dimentions: (1 x 1 x wavenumber)
%   -alpha_O2:[1/m] (range x time) calculated absorption from data
%   -SNRm:[none] (range x time) calculated SNR mask
%   -cloud_SDm_above:[none] (range x time) calculated cloud mask
%
%Outputs:
%   -T_final: [K] (range x time x iteration) Final calculated temperature profile
%   -Lapse: [K/m] (range x time x iteration) Fitted lapse rate
%   -Ts_fit: [K] (range x time x iteration) Fitted surface temperature
%   -P_final: [atm] (range x time x iteration) Final calculated pressure
%   -mean_lapse_rate: [atm] (range x time x iteration) Mean lapse rate for
%   whole day
%   -exclusion: [none] (range x time x iteration) Points excluded from lapse rate fit


%==Constants
g0 = 9.80665;               %[m/s/s]Gravitational acceleration
M_air = 0.0289644;          %[kg/mol] Molar mass of air
q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio 
kB = 1.38065E-23;           %[J/K] Boltzman's constant 
%No = 2.47937E25;            % Loschmidt's number [1/m^3] (referenced to 296 K and 1 atm)
N_A = 6.02214076e23;        %[1/mol] Avagadro's number
R = kB * N_A;               %[J/K/mol] universal gas constant
h = 6.626E-34;              %[J s] Planck's constant
c = 2.99792458E8;           %[m/s] speed of light 

%==O2 line parameters
S0_O2 = 4.889e-26;          %[cm-1/molec cm-2] line strength
S0_O2 = S0_O2 / 100;        %[m molecule-1] absoption line strength at T=296K
E_lower = 1420.7631;        %[cm-1] ground state energy
E_lower = E_lower * 100;    %[m-1] ground state energy
ep = E_lower*h*c;           %[J]

T0 = 296;                   %[K] reference temperature

loop = 60;%number of times to do iterative temperature retrieval loop

del_r = rm(2)-rm(1);        %[m]set range difference
Ts = T(1,:);                %[K] surface temp
Ps = P(1,:);                %[atm] surface press

gamma = g0 * M_air / R;     %[K/m]gravity molar mass of air and gas constant

%Set inital temperature guess
Tg = T;                     %[K]

%Preallocate memory
deltaT = zeros(length(rm),length(ts),loop);
T_ret = zeros(length(rm),length(ts),loop);
exponent = zeros(length(rm),length(ts),loop);
P_ii = zeros(length(rm),length(ts));
Lapse = zeros(length(rm),length(ts),loop);
Ts_fit = zeros(length(rm),length(ts),loop);
exclusion = zeros(length(rm),length(ts),loop);

logicalExc = true(length(rm),length(ts));

starting_lapse_rate = -0.0065;%[K/m] typical lapse rate

mean_lapse_rate = ones(1,1,1);
mean_lapse_rate(:,:,:) = starting_lapse_rate;%[K/m] variable to set to for unfitted points

L_t_window = 2;

%Fit lapse rate bounds
lower_alt_threshold = 1000;  %[m] only preform a fit above
upper_lapse_bound = -0.004;  %[K/m]
lower_lapse_bound = -0.010;  %[K/m]

fprintf('Loop')
for i = 1:loop
    fprintf(' %d',i)
    %calculate lapse rate
    for time = 1:length(ts)
        % Set fit exclusion zones
        exclusion(:,time,i) = rm<lower_alt_threshold | SNRm(:,time)==0 | cloud_SDm_above(:,time) ~= 1;
        % Do not fit if the temperature guess is all nans
        if isnan(Tg(:,time))
            % Set fits to mean lapse rates
            Lapse(:,time,i) = mean_lapse_rate(:,:,1);%set lapserate to default
            Ts_fit(:,time,i) = Ts(time);   
            
        %make sure there are enough points to fit
        elseif sum(1-exclusion(:,time,i)) > 20 
            %fo = fitoptions('Exclude',exclusion(:,time,i),...
                            %'StartPoint',[starting_lapse_rate,Ts(time)]);%fit options
% % % %             Lapse_fit = fit(rm,Tg(:,time),polyf,...
% % % %                             'Exclude',exclusion(:,time,i),...
% % % %                             'Lower',[lower_lapse_bound,lower_Ts_bound(time)],...
% % % %                             'Upper',[upper_lapse_bound,upper_Ts_bound(time)],...
% % % %                             'Weights',o2on_SNR(:,time));
% % % %                                         
% % % %             Lapse(:,time,i) = Lapse_fit.p1;
% % % %             Ts_fit(:,time,i) = Lapse_fit.p2;
            
            %Custom least squares for speed
            %initalize logical and V
            if i==1   
                logicalExc(:,time) = ~logical(exclusion(:,time,i)); 
            end
            V = [rm(logicalExc(:,time)) ones(length(rm(logicalExc(:,time))),1)];
            %calculate polynomial coefficients     
            p = V\Tg(logicalExc(:,time),time);
            
            Lapse(:,time,i) = p(1);
            Ts_fit(:,time,i) = p(2);
            
            % Keep lapse rate within limits
            if Lapse(:,time,i) <lower_lapse_bound
                Lapse(:,time,i) = lower_lapse_bound;
            elseif Lapse(:,time,i) >upper_lapse_bound
                Lapse(:,time,i) = upper_lapse_bound;
            end
            
        %if not enough points, set the first datapoints to the mean lapse rate 
        elseif time <= L_t_window/2
            disp('time <=L_t_window')
            Lapse(:,time,i) = mean_lapse_rate(:,:,1);%set lapserate to default
            Ts_fit(:,time,i) = Ts(time);
        %if not enough points set lapse to mean of lastfew means
        else 
            Lapse(:,time,i) = mean(Lapse(:,time-L_t_window/2:time-1,i),2);
            Ts_fit(:,time,i) = Ts(time);
            
         end
    end
    
    %calclate total mean
    mean_lapse_rate(:,:,1) = mean(Lapse(1,:,i),2);    
    
    % === Pressure Profile ===
    P_i = Ps;
    T_i = Tg;
    L_i = Lapse(:,:,i);
    
    for j = 1:length(rm)
        P_ii(j,:) = P_i.*(T_i(j,:)./(T_i(j,:) + L_i(j,:).*del_r)).^(gamma./L_i(j,:));
        P_i = P_ii(j,:);
    end
    
    Pg = [Ps; P_ii(1:end-1,:)];%concatenate arrays to set ground pressure to measured value
    
    %calculate air number density
    n_L = Pg .* 101325 / kB ./ Tg;          %[1/m^3] Loschmidt's number; with pressure converted to Pa
    q_WV = WV ./ n_L;                       %[unitless] mixing ratio of water vapor based from water vapor number density 
    q = q_O2 .* (1 - q_WV);
    
    %update lineshape function
    [~,~,g,Line] = absorption_O2_770_model(Tg,Pg,nu_scan,WV);%[m] lineshape function 
    
    %Calculate Coefficients
    exponent(:,:,i) = gamma ./ Lapse(:,:,i);
    C1 = S0_O2 * T0 * (Pg*101325) * exp(ep/kB/T0) ./ (kB * Tg .^(-exponent(:,:,i)));%[K/(molec*m^2)]pressure converted to Pa
    C2 = Tg .^ (-exponent(:,:,i) - 2) .* exp(-ep/kB./Tg);%[K^-1]
    C3 = (-exponent(:,:,i) - 2) ./ Tg + ep./(kB.*Tg.^2);%[K^-1]
    %Calculate change in temperature from last
    %deltaT(:,:,i) = (alpha_O2 - C1 .* C2 .* g .* q) ./ (C1 .* C2 .* C3 .* g .* q); %[K] calculate a change in temperatre
    deltaT(:,:,i) = (alpha_O2 - C1 .* C2 .* Line{2}.lineshape .* q) ./ (C1 .* C2 .* C3 .* Line{2}.lineshape .* q); %[K] calculate a change in temperatre

    %Line coefficients
    epa = Line{1}.E_lower*h*c;
    C1a = Line{1}.S0 * T0 * (Pg*101325) * exp(epa/kB/T0) ./ (kB * Tg .^(-exponent(:,:,i)));%[K/(molec*m^2)]pressure converted to Pa
    C2a = Tg .^ (-exponent(:,:,i) - 2) .* exp(-epa/kB./Tg);%[K^-1]
    C3a = (-exponent(:,:,i) - 2) ./ Tg + epa./(kB.*Tg.^2);%[K^-1]
    epb = Line{2}.E_lower*h*c;
    C1b = Line{2}.S0 * T0 * (Pg*101325) * exp(epb/kB/T0) ./ (kB * Tg .^(-exponent(:,:,i)));%[K/(molec*m^2)]pressure converted to Pa
    C2b = Tg .^ (-exponent(:,:,i) - 2) .* exp(-epb/kB./Tg);%[K^-1]
    C3b = (-exponent(:,:,i) - 2) ./ Tg + epb./(kB.*Tg.^2);%[K^-1]
    %Calculate change in temperature from last
    deltaT(:,:,i) = (alpha_O2 - C1a.*C2a.*Line{1}.lineshape.*q - C1b.*C2b.*Line{2}.lineshape.*q)./(C1a.*C2a.*C3a.*Line{1}.lineshape.*q+C1b.*C2b.*C3b.*Line{2}.lineshape.*q); %[K] calculate a change in temperatre

    % Limit deltaT to plus or minus 2 K
    deltaT(deltaT > 2) = 2;
    deltaT(deltaT < -2) = -2;

    % Update temperature profile guress
    T_ret(:,:,i) = Tg + deltaT(:,:,i);            
    Tg = fillmissing(T_ret(:,:,i),'nearest',1);
    Tg = fillmissing(Tg,'nearest',2);
    
    %Fill data points above cloud and snr mask with lapse rate values
    tempT = Ts_fit(:,:,i) + rm.* Lapse(:,:,i);
    Tg(SNRm(:,:)==0 | cloud_SDm_above(:,:) ~= 1) = tempT(SNRm(:,:)==0 | cloud_SDm_above(:,:) ~= 1);
    
end

fprintf('\n')
% Final temperature and pressure
T_final = Tg;
P_final = Ps .* (Ts./(Ts+Lapse(:,:,end).*(rm-60))).^(0.0341623./Lapse(:,:,end));

exclusion = exclusion(:,:,end);
end