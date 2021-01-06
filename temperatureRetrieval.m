function [T_final,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion] =  temperatureRetrieval(T,ts,rm,P,WV,nu_scan,alpha_O2,SNRm,cloud_SDm_above,o2on_SNR)
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
%   -o2on_SNR:[none] (range x time) signal to noise ratio for online counts
%
%Outputs:
%   -T_final: [K] (range x time x iteration) Final calculated temperature profile
%   -Lapse: [K/m] (range x time x iteration) Fitted lapse rate
%   -Ts_fit: [K] (range x time x iteration) Fitted surface temperature
%   -P_final: [atm] (range x time x iteration) Final calculated pressure
%   -mean_lapse_rate: [atm] (range x time x iteration) Mean lapse rate for
%   whole day
%   -exclusion: [none] (range x time x iteration) Points excluded from lapse rate fit


%constants
g0 = 9.80665;               %[m/s/s]Gravitational acceleration
M_air = 0.0289644;          %[kg/mol] Molar mass of air
q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio 
kB = 1.38065E-23;           %[J/K] Boltzman's constant 
%No = 2.47937E25;            % Loschmidt's number [1/m^3] (referenced to 296 K and 1 atm)
N_A = 6.02214076e23;        %[1/mol] Avagadro's number
R = kB * N_A;               %[J/K/mol] universal gas constant


%sw_single = 1.106e-25;%[cm-1/molec cm-2] line strength
%S0_O2 = 1.106e-25;          %[cm molecule-1] line strength
S0_O2 = 4.889e-26;
S0_O2 = S0_O2 / 100;        %[m molecule-1] absoption line strength at T=296K
%S0_O2 = S0_O2 * 0.75;
%E_lower = 1248.2013;        %[cm-1] ground state energy
E_lower = 1420.7631;
E_lower = E_lower * 100;    %[m-1] ground state energy
h = 6.626E-34;              %[J s] Planck's constant
c = 2.99792458E8;           %[m/s] speed of light 
ep = E_lower*h*c;           %[J]

T0 = 296;%[K] reference temperature

L_t_window = 30; %[min] time moving average
L_r_window = 8; %[x37.5m] range moving average
L_k = ones(L_r_window,L_t_window)./(L_r_window*L_t_window); %lapse rate kernel


loop = 20;%number of times to do iterative temperature retrieval loop

del_r = rm(2)-rm(1);%[m]set range difference
t_step = ts(2)-ts(1);%[s] time step
Ts = T(1,:);%[K]
Ps = P(1,:);%[atm]

gamma = g0 * M_air / R; %[K/m]gravity molar mass of air and gas constant

%Set inital temperature guess
Tg = T;%[K]

%Preallocate memory
L_fit = zeros(length(rm),length(ts),loop);
L_fit_sm = zeros(length(rm),length(ts),loop);
L_mean = zeros(1,length(ts),loop);
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

%fit lapse rate
lower_alt_threshold = 1100;%[m] only preform a fit above
lower_alt_threshold = 1500;%[m] only preform a fit above
lower_alt_threshold = 1000;%[m] only preform a fit above

%lapse rate bounds
upper_lapse_bound = -0.004; %[K/m]
lower_lapse_bound = -0.010;  %[K/m]
%lower_lapse_bound = -0.05;  %[K/m]


upper_Ts_bound = Ts + 20;%[K]
lower_Ts_bound = Ts - 20;%[K]

%Set fittype, reduces computation time by a lot outside loop
polyf = fittype('poly1');

fprintf('Loop')
for i = 1:loop
    fprintf(' %d',i)
    %calculate lapse rate
    %[~,L] = gradient(Tg,ts,rm);             %[K/m] Lapse rate from temperature guess with spacing of range vector
     
    for time = 1:length(ts)
        exclusion(:,time,i) = rm<lower_alt_threshold | SNRm(:,time)==0 | cloud_SDm_above(:,time) ~= 1;
        if isnan(Tg(:,time))
            %Tg(1,time);
            %time;
            Lapse(:,time,i) = mean_lapse_rate(:,:,1);%set lapserate to default
            Ts_fit(:,time,i) = Ts(time);
        elseif sum(1-exclusion(:,time,i)) > 20 & isnan(o2on_SNR(1,time)) %make sure there are enough points to fit
            %fo = fitoptions('Exclude',exclusion(:,time,i),...
                            %'StartPoint',[starting_lapse_rate,Ts(time)]);%fit options
                            %Ts(time);
                            %[lower_lapse_bound,lower_Ts_bound(time)];
            Lapse_fit = fit(rm,Tg(:,time),polyf,...
                            'Exclude',exclusion(:,time,i),...
                            'Lower',[lower_lapse_bound,-inf],...
                             'Upper',[upper_lapse_bound,inf]);
%                             'Lower',[lower_lapse_bound,lower_Ts_bound(time)],...
%                             'Upper',[upper_lapse_bound,upper_Ts_bound(time)]);


                                        
            Lapse(:,time,i) = Lapse_fit.p1;
            Ts_fit(:,time,i) = Lapse_fit.p2;
            
        elseif sum(1-exclusion(:,time,i)) > 20 %make sure there are enough points to fit
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
            
            %Custom least squares
            
            %V = [Tg(:,time) ones(length(Tg(logicalExc,time)),1)];
% %             V = [Tg(logicalExc,time) ones(length(Tg(logicalExc,time)),1)];
% %             p = V\rm(logicalExc);

            %initalize logical and V
            if i==1
                
                logicalExc(:,time) = ~logical(exclusion(:,time,i));
                %V = ones(sum(logicalExc(:,time)),2,length(ts));
                
            end
            V = [rm(logicalExc(:,time)) ones(length(rm(logicalExc(:,time))),1)];
            %calculate polynomial coefficients
            %V(:,:,time)
            
            p = V\Tg(logicalExc(:,time),time);
            
            Lapse(:,time,i) = p(1);
            Ts_fit(:,time,i) = p(2);
            
            %p2 = polyfit(rm(logicalExc),
            
           
        elseif time == 1
            Lapse(:,time,i) = mean_lapse_rate(:,:,1);%set lapserate to default
            Ts_fit(:,time,i) = Ts(time);
        elseif time <= L_t_window/2
            Lapse(:,time,i) = mean_lapse_rate(:,:,1);%set lapserate to default
            Ts_fit(:,time,i) = Ts(time);
        else 
            Lapse(:,time,i) = mean(Lapse(:,time-L_t_window/2:time-1,i),2);
            Ts_fit(:,time,i) = Ts(time);
            
%         else
%             Lapse(:,time,i) = mean_lapse_rate(:,:,i);%set lapserate to previous lapse rate
%             Ts_fit(:,time,i) = Ts(time);
         end
    end
    
    mean_lapse_rate(:,:,1) = mean(Lapse(1,:,i),2);
    
    % Smooth over range and time
   % L_fit(:,:,i) = L;
    %L_fit_offset = diff(Tg)./del_r;
    %L_fit(:,:,i) = interp1(rm(2:end) - (del_r/2),L_fit_offset,rm,'nearest','extrap');
    %L_fit_sm(:,:,lp-1) = movmean(L_fit(:,:,lp-1),L_t_window,2);
   % L_fit_pad = padarray(L_fit(:,:,i),[L_r_window/2,L_t_window/2],'replicate');
   % L_fit_filt = filter2(L_k,L_fit_pad,'valid');
    %L_fit_intp = interp2(ts-t_step/2,rm-del_r/2,L_fit_filt(1:end-1,1:end-1),ts,rm);
   % L_fit_fill = fillmissing(L_fit_intp,'nearest',1);        % Fill in NaNs in dimension 1
    %L_fit_sm(:,:,i) = fillmissing(L_fit_fill,'nearest',2);   % Fill in NaNs in dimension 2

    % Mean lapse rate (average over range)
    %L_mean(:,:,i) = mean(L_fit_sm(:,:,i) .* (1 - exclusion(:,:,i)));
    
    % Use this to get back to temperature from lapse rate:
    % tg_fit = T_surf + cumtrapz((rm-60),L_fit,1);
    
    
    % === Pressure Profile ===
    P_i = Ps;
    T_i = Tg;
    %L_i = L_fit(:,:,i);
    L_i = Lapse(:,:,i);
    
    for j = 1:length(rm)
    P_ii(j,:) = P_i.*(T_i(j,:)./(T_i(j,:) + L_i(j,:).*del_r)).^(gamma./L_i(j,:));
    P_i = P_ii(j,:);
    end
    
    Pg = [Ps; P_ii(1:end-1,:)];%concatenate arrays to set ground pressure to measured value
    
    %Pressure
    %Pg = Ps .* (Ts ./ Tg) .^ (gamma ./ Lapse(:,:,i));  %[atm]
   
    
    %calculate air numer density
    n_L = Pg .* 101325 / kB ./ Tg;          %[1/m^3] Loschmidt's number; with pressure converted to Pa
    q_WV = WV ./ n_L;                       %[unitless] mixing ratio of water vapor based from water vapor number density
    n_O2 = q_O2 .* (1 - q_WV) .* n_L;       %[1/m^3] number density of O2
    
    q = q_O2 .* (1 - q_WV);
    
    %update lineshape function
    [~,~,f] = absorption_O2_770_model(Tg,Pg,nu_scan,WV);
    %[~,~,f] = absorption_O2_770_model_wavenumber(Tg,Pg,nu_scan,WV);
    g = f;                                  %[m] lineshape function 
    
    %exponent(:,:,i) = gamma ./ L_fit(:,:,i);
    exponent(:,:,i) = gamma ./ Lapse(:,:,i);
    C1 = S0_O2 * T0 * (Pg*101325) * exp(ep/kB/T0) ./ (kB * Tg .^(-exponent(:,:,i)));%[K/(molec*m^2)]pressure converted to Pa

    C2 = Tg .^ (-exponent(:,:,i) - 2) .* exp(-ep/kB./Tg);%[K^-1]

    C3 = (-exponent(:,:,i) - 2) ./ Tg + ep./(kB.*Tg.^2);%[K^-1]

    deltaT(:,:,i) = (alpha_O2 - C1 .* C2 .* g .* q) ./ (C1 .* C2 .* C3 .* g .* q); %[K] calculate a change in temperatre

    % Limit deltaT to plus or minus 2 K
    deltaT(deltaT > 2) = 2;
    deltaT(deltaT < -2) = -2;

    %update temperature profile
    T_ret(:,:,i) = Tg + deltaT(:,:,i);            
    Tg = fillmissing(T_ret(:,:,i),'nearest',1);
    Tg = fillmissing(Tg,'nearest',2);
    
    %Fill data points above cloud and snr mask with lapse rate values
    tempT = Ts_fit(:,:,i) + rm.* Lapse(:,:,i);
    Tg(SNRm(:,:)==0 | cloud_SDm_above(:,:) ~= 1) = tempT(SNRm(:,:)==0 | cloud_SDm_above(:,:) ~= 1);
    

end

fprintf('\n')

new_T_pad = padarray(Tg,[L_r_window/2,L_t_window/2],'replicate');
new_T_filt = filter2(L_k,new_T_pad,'valid');
new_T_smooth = interp2(ts-L_t_window/2,rm-L_r_window/2,new_T_filt(1:end-1,1:end-1),ts,rm);
new_T_smooth = fillmissing(new_T_smooth,'nearest',1); % Fill in NaNs in dimension 1
new_T_smooth = fillmissing(new_T_smooth,'nearest',2); % Fill in NaNs in dimension 2

T_final = new_T_smooth;
T_final_test = T_final;
%T_final = Tg;
P_final = Ps .* (Ts./(Ts+Lapse(:,:,end).*(rm-60))).^(0.0341623./Lapse(:,:,end));

exclusion = exclusion(:,:,end);
end