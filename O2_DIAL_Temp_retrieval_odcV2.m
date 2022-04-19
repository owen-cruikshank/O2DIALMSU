% O2_DIAL_Temp_retrieval_odcV2.m
% Analysis program for O2 DIAL instrument from 
clear 

%=Date Range

date_start = datetime(2021,10,22,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2021,10,23,'TimeZone','UTC');%yyyy,mm,dd

span_days = date_start:date_end;

%=Time and range averaging
Options.intTime = 5;  %[min] Integration time
Options.intRange = 1; %[bins] Integration range

Options.t_avg = 1;     %[bins] Time smoothing bins
Options.oversample = 1; %[bins] Range smoothing bins

%====================
%==== Constants =====
%====================
Constant.g0 = 9.80665;       %[m/s^2] Gravitational acceleration 
Constant.M_air = 0.0289644;  %[kg/mol] Molar mass of Earth's air 
Constant.R = 8.3144598;      %[J/(mol*K)] Universal gas constant 

Constant.c = 2.99792458E8;           %[m/s] Speed of light 
Constant.kb = 1.38065E-23;           %[J/K][m^2 kg s-2 K-1] Boltzman's constant 
Constant.h = 6.626E-34;              %[Js] Planck's constant 
Constant.mo2 = 5.314E-26;            %[kg] Mass O2 molecule 
Constant.mWV = 2.9915e-26;           %[kg] Mass H2O molecule
Constant.m_air = 4.792E-26;          %[kg] Mass of air
Constant.q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
Constant.No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)

%======================
%==== Load data
%==== photon counts, radiosondes, weather station data, HSRL data, etc.
%======================

%==Load data from MSU instument==
%[Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadMSUdata(span_days,Options,Constant);
%   HSRL.BSR = ones(size(Counts.o2on));
%   HSRL.Ba = zeros(size(HSRL.Ba));
%  HSRL.BSR(isnan(Counts.o2on))=nan;

% ==Load data from SGP ARM instument==
%date_begin = datetime(2019,4,17); date_end   = datetime(2019,4,22);
%span_days = date_begin:date_end;        % Days in set [datetime]
%[Range,Time,Counts,Sonde,Model,Spectrum,BSR,Data] = loadSGPdata(span_days,Options,Constant);

%==Load data from Boulder instument==
% date_begin = datetime(2020,8,29); 
% date_end   = datetime(2020,9,8);

date_begin   = datetime(2021,6,20);
date_end   = datetime(2021,6,22);
span_days = date_begin:date_end;        % Days in set [datetime]

Options.MPDname = '05';
%Options.MPDname = '01';
[Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] =loadBoulderdata3(span_days,Options,Constant);

%%
%======================
%= Cloud and SNR mask =
%======================
disp('Calculating masks')
cloud_p_point = 113.9183; %(hours) Time to plot mask data
SNR_threshold = 50;
SD_threshold = 5;
BGmult =1; % Multiplier for background for SNR calculation
lowAlt = 400;
%[SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2_counts(Counts.o2on,Counts.o2off,Range.rm,Time.ts,cloud_p_point,SNR_threshold,SD_threshold,Options.oversample,Options.t_avg,Counts,BGmult);
[SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2_BSR(cloud_p_point,SNR_threshold,SD_threshold,Options.oversample,Options.t_avg,Counts,BGmult,Time,Range,HSRL.BSR,lowAlt);

%%
% ==== Afterpulse correction ====
% pulseON = 4.146e+07*10*(Range.rm).^(-2.536*1);
% % pulseOFF = 4.146e+07*1*(rm_over).^(-2.536*1.0);
% pulseOFF = 3.682e+07*10*1.3*(Range.rm).^(-2.512*1);
% 
% % pulseON = 4.146e+07*10*1.3*(Range.rm).^(-2.536*1);
% % % pulseOFF = 4.146e+07*1*(rm_over).^(-2.536*1.0);
% % pulseOFF = 3.682e+07*10*1*(Range.rm).^(-2.512*1);
% Counts.o2on = Counts.o2on-pulseON;
% Counts.o2off = Counts.o2off-pulseOFF;

%%
%%%for aaa = 1:2
%     if aaa >1
%         Tnew = Temperature.T_finalm;
%         Pnew = Temperature.Patm_final;
%         Tnew(isnan(Tnew)) = Model.T(isnan(Tnew));
%         Model.T = Tnew;
%     end

%%
%==========================
%= Pertabative absorption =
%==========================
disp('Calculating absorption')
% ==== Afterpulse correction for absorption
%for ii=1
    %if ii==2
    %    ln_o2 = log(((Counts.o2on(ind_r_lo,:)-delta(ind_r_lo,:)).*(Counts.o2off(ind_r_hi,:)-delta(ind_r_lo,:)))./((Counts.o2on(ind_r_hi,:)-delta(ind_r_lo,:)).*(Counts.o2off(ind_r_lo,:)+delta(ind_r_lo,:)))); % Natural log of counts
    %else
%         ind_r_lo = 1:Range.i_range-Options.oversample;                                            % High range vector
%         ind_r_hi = 1+Options.oversample:Range.i_range;                                            % Low range vector
%         ln_o2 = log((Counts.o2on(ind_r_lo,:).*Counts.o2off(ind_r_hi,:))./(Counts.o2on(ind_r_hi,:).*Counts.o2off(ind_r_lo,:))); % Natural log of counts
%     
        %load('Delta.mat')
%         load('Delta2.mat')
%  load('Delta3.mat')
%         % Delta=Delta2;
%          Delta=Delta(:,611);
%          binShift =5;
%          Delta = [Delta(binShift:end,:); zeros(length(Range.rm)-size(Delta(binShift:end,:),1),1)];
% load('AfterPulse.mat')
% Delta = pulseON(4:133+3);
% % Delta = pulseON(1:133);
%           ln_o2 = log(((Counts.o2on(ind_r_lo,:)-Delta(ind_r_lo,:)).*(Counts.o2off(ind_r_hi,:)-Delta(ind_r_hi,:)))./((Counts.o2on(ind_r_hi,:)-Delta(ind_r_hi,:)).*(Counts.o2off(ind_r_lo,:)-Delta(ind_r_lo,:)))); % Natural log of counts
        %end
        
% === Zeroth Order ===
ind_r_lo = 1:Range.i_range-Options.oversample;                                            % High range vector
ind_r_hi = 1+Options.oversample:Range.i_range;                                            % Low range vector
ln_o2 = log((Counts.o2on(ind_r_lo,:).*Counts.o2off(ind_r_hi,:))./(Counts.o2on(ind_r_hi,:).*Counts.o2off(ind_r_lo,:))); % Natural log of counts
    
Alpha.alpha_0_raw = ln_o2./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
Alpha.alpha_0 = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0_raw,Time.ts,Range.rm);

Alpha.alpha_0 = Alpha.alpha_0 - Model.absorption_off; %Subtract offline absorption

% fon = Counts.fon-Counts.fon_bg;
% foff = Counts.foff-Counts.foff_bg;
% gon = Counts.gon-Counts.fon_bg;
% goff = Counts.goff-Counts.foff_bg;
% ln_o2 = log((fon(ind_r_lo,:).*foff(ind_r_hi,:))./(fon(ind_r_hi,:).*foff(ind_r_lo,:))); % Natural log of counts  
% Alpha.alpha_0f_raw = ln_o2./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
% Alpha.alpha_0f = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0f_raw,Time.ts,Range.rm);
% Alpha.alpha_0f = real(Alpha.alpha_0f);
% 
% ln_o2 = log((gon(ind_r_lo,:).*goff(ind_r_hi,:))./(gon(ind_r_hi,:).*goff(ind_r_lo,:))); % Natural log of counts  
% Alpha.alpha_0g_raw = ln_o2./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
% Alpha.alpha_0g = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0g_raw,Time.ts,Range.rm);
% Alpha.alpha_0g = real(Alpha.alpha_0g);

ln_o2 = log((Counts.wvon(ind_r_lo,:).*Counts.wvoff(ind_r_hi,:))./(Counts.wvon(ind_r_hi,:).*Counts.wvoff(ind_r_lo,:))); % Natural log of counts  
Alpha.alpha_0wv_raw = ln_o2./2./(Range.rangeBin*Options.oversample); 
Alpha.alpha_0wv = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0wv_raw,Time.ts,Range.rm);

clear ln_o2 foff fon goff gon

[~,cross_section,~,~] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvon,Alpha.alpha_0wv);

[~,cross_sectionOff,~,~] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvoff,Alpha.alpha_0wv);

N_wv0 = Alpha.alpha_0wv./(cross_section-cross_sectionOff);
%%
%%==============Set Model WV to DIAL======

%%%%Purturbative WV
load('TransmittanceData.mat')
T_etalon_on = double(interp1(double(OnlineWavelength)*10^9,OnlineCombinedTransmittance,Spectrum.lambda_scan_3D_short));
T_etalon_off = double(interp1(double(OfflineWavelength)*10^9,OfflineCombinedTransmittance,Spectrum.lambda_scan_3D_short_off));

altitude = 1.5;%altitude in km
altitude = 1.5719;

[Alpha.alpha_1wv, Alpha.alpha_2wv,Spectrum] = pertAbsorptionwv(Alpha.alpha_0wv, T_etalon_on, Model, Range, Time, Spectrum, HSRL, ind_r_lo,ind_r_hi, Options, Constant, altitude);

N_wv = (Alpha.alpha_0wv+Alpha.alpha_1wv+ Alpha.alpha_2wv)./(cross_section-cross_sectionOff);

%Smooth WV
k = ones(4,8)./(4*8);     % Kernel
%k = ones(3,7)./(3*7);     % Kernel



N_wvm = nanconv(N_wv,k,'edge','nanout');
N_wvm(cloud_SDm_above)=nan;

N_wv0m = nanconv(N_wv0,k,'edge','nanout');
N_wv0m(cloud_SDm_above)=nan;

AbsHumm = N_wvm.*Constant.mWV*1000; %[g/m3]
AbsHum0m = N_wv0m.*Constant.mWV*1000; %[g/m3]

AbsHumRawm = N_wv.*Constant.mWV*1000; %[g/m3]
AbsHumRawm(cloud_SDm_above)=nan;
AbsHum0Rawm = N_wv0.*Constant.mWV*1000; %[g/m3]
AbsHum0Rawm(cloud_SDm_above)=nan;
%%%%%%SET MODEL TO WV RETRIEVAL
%%%%%%%%%%%Model.WV = fillmissing(N_wvm,'linear');


%%
% === Total alpha ===

Alpha.alpha_total_rawwv = Alpha.alpha_0wv + Alpha.alpha_1wv + Alpha.alpha_2wv;

N_wv = Alpha.alpha_total_rawwv./cross_section; %[molecule/m3] wv number density



% === Smoothing zeroth order

Alpha.alpha_0=real(Alpha.alpha_0);



% === Molecular alpha calcultion =====


%%
% === Purtabative absorption ===
%loading etalon transmission

[Alpha.alpha_total_raw,Alpha.alpha_1,Alpha.alpha_2,Spectrum] = pertAbsorption(Alpha.alpha_0, T_etalon_on, Model, Range, Time, Spectrum, HSRL, ind_r_lo,ind_r_hi, Options);

%%%[Alpha.alpha_total_rawf,Alpha.alpha_1f,Alpha.alpha_2f,~] = pertAbsorption(Alpha.alpha_0f, T_etalon_on, Model, Range, Time,Spectrum, HSRL, ind_r_lo,ind_r_hi, Options);
%%%[Alpha.alpha_total_rawg,Alpha.alpha_1g,Alpha.alpha_2g,~] = pertAbsorption(Alpha.alpha_0g, T_etalon_on, Model, Range, Time,Spectrum, HSRL, ind_r_lo,ind_r_hi, Options);
%%



%== Force total alpha to its modeled surface value ==
%[~,cut] = min(abs(Range.rm-500));             % Index where rm is closest to chosen value
cut=2;
Alpha.alpha_total_cut = [Model.absorption(1,:); NaN((cut - 2),Time.i_time); Alpha.alpha_total_raw(cut:end,:)];
Alpha.alpha_total_cut = fillmissing(Alpha.alpha_total_cut,'linear');

%%
%===== Soothing alpha =====
%%%k = ones(10,4)./(10*4);     % Kernel
% k = ones(4,6)./(4*6);     % Kernel
% k = ones(8,6)./(8*6);     % Kernel
% k = ones(8,8)./(8*8);     % Kernel
% k = ones(2,4)./(2*4);
k = ones(1,1)./(1*1);

gg = cloud_SDm_above.*SNRm;
gg(gg<=0)=nan;
%Appy cloud mask before smoothing
Alpha.alpha_total_cut(isnan(gg)) = NaN;          % Replace mask with NaNs

Alpha.alpha_total_rawmf = Alpha.alpha_total_rawf;
Alpha.alpha_total_rawmg = Alpha.alpha_total_rawg;
Alpha.alpha_total_rawmf(isnan(gg))=nan;
Alpha.alpha_total_rawmg(isnan(gg))=nan;

%%%%[Ez,Et,minSigz,minSigt] = findMinE(Alpha.alpha_total_rawmf,Alpha.alpha_total_rawmg,0);
%%%%[Alpha.alpha_total_filt2] = applyFilter(minSigz,minSigt,Alpha.alpha_total_cut);

%%

%Convolve with kernel
Alpha.alpha_total_filt = nanconv(Alpha.alpha_total_cut,k,'edge','nanout');
Alpha.alpha_totals = Alpha.alpha_total_filt;

%%%Alpha.alpha_totals = Alpha.alpha_total_filt2;

Alpha.alpha_total_filtf = nanconv(Alpha.alpha_total_rawf,k,'edge','nanout');
Alpha.alpha_totalsf = Alpha.alpha_total_filtf;
Alpha.alpha_total_filtg = nanconv(Alpha.alpha_total_rawg,k,'edge','nanout');
Alpha.alpha_totalsg = Alpha.alpha_total_filtg;


%%
% apply SNR mask again
Alpha.alpha_0m = Alpha.alpha_0;
Alpha.alpha_0m(isnan(gg)) = NaN;                  % Replace mask with NaNs

Alpha.alpha_total = real(Alpha.alpha_totals);
Alpha.alpha_totalm = Alpha.alpha_total;
Alpha.alpha_totalm(isnan(gg)) = NaN;          % Replace mask with NaNs


%Smooth WV
gg = cloud_SDm_above.*SNRm;
gg(gg<=0)=nan;
N_wvm = nanconv(N_wv,k,'edge','nanout');
N_wvm(isnan(gg))=nan;

N_wv0m = nanconv(N_wv0,k,'edge','nanout');
N_wv0m(isnan(gg))=nan;

AbsHumm = N_wvm.*Constant.mWV*1000; %[g/m3]
AbsHum0m = N_wv0m.*Constant.mWV*1000; %[g/m3]

AbsHumRawm = N_wv.*Constant.mWV*1000; %[g/m3]
AbsHumRawm(isnan(gg))=nan;
AbsHum0Rawm = N_wv0.*Constant.mWV*1000; %[g/m3]
AbsHum0Rawm(isnan(gg))=nan;



%%
%===== Find afterpulsing correction ======
%delta = findDelta(Model.absorption(:,:),alpha_total(:,:),Counts.o2on(:,:),Counts.o2off(:,:),Range.rm,Time.ts);
Counts.delta = findDelta(Model.absorption(:,:),Alpha.alpha_total_raw(:,:),Counts.o2on(:,:),Counts.o2off(:,:),Range.rm,Time.ts);
Counts.delta=-Counts.delta;

Counts.delta = [Counts.delta; zeros(length(Range.rm)-size(Counts.delta,1),length(Time.ts))];
Counts.delta = Counts.delta/4;

% % deltaModel = findDelta(Model.absorption(:,:),alpha_total_raw(:,:),Model.N_on(:,:),Model.N_off(:,:),Range.rm,Time.ts);
% % deltaModel=-deltaModel;
% % deltaModel = [deltaModel; zeros(length(Range.rm)-size(deltaModel,1),length(Time.ts))];
% % 
% % % binShift = 6;
% % % delta = [delta(binShift:end,:); zeros(length(Range.rm)-size(delta(binShift:end,:),1),length(Time.ts))];
% % 
% % alpha_d = 1/2/(Range.rangeBin).*log((Counts.o2off(2:end,:)-delta(2:end,:)).*(Counts.o2on(1:end-1,:)-delta(1:end-1,:))./(Counts.o2off(1:end-1,:)-delta(1:end-1,:))./(Counts.o2on(2:end,:)-delta(2:end,:)));
% % 
% % alpha_dModel = 1/2/(Range.rangeBin).*log((Counts.o2off(2:end,:)-deltaModel(2:end,:)).*(Counts.o2on(1:end-1,:)-deltaModel(1:end-1,:))./(Counts.o2off(1:end-1,:)-deltaModel(1:end-1,:))./(Counts.o2on(2:end,:)-deltaModel(2:end,:)));
% % 
% % %alpha_d_pad = padarray(alpha_d,[8/2,6/2],'replicate');
% % %alpha_d_filt = filter2(k,alpha_d_pad,'valid');
% % %alpha_d = interp2(Time.ts-Time.t_step/2,Range.rm(1:end-1)-Range.rangeBin/2,alpha_d_filt(1:end-1,1:end-1),Time.ts,Range.rm);
% % alpha_d_filt = filter2(k,alpha_d,'same');

%%
%==========================
%= Temperature Function =
%==========================
disp('Temp retrieval function')

% Testing different models
%Model.T = Model.Ts + lapseRate .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
%  Model.T = 240*ones(Range.i_range,Time.i_time);
%  Model.Ts = Model.T(1,:);
% Model.P = Model.Ps .* (Model.Ts./Model.T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  

tic
[Temperature.T_final_test,Temperature.L_fit_sm_test,Temperature.Ts_fit,Temperature.Patm_final,Temperature.mean_lapse_rate,Temperature.exclusion,Temperature.Titer] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_totals,SNRm,cloud_SDm_above);

[Temperature.T_final_testf,Temperature.L_fit_sm_test,Temperature.Ts_fit,Temperature.Patm_final,Temperature.mean_lapse_rate,Temperature.exclusion,Temperature.Titer] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_totalsf,SNRm,cloud_SDm_above);
[Temperature.T_final_testg,Temperature.L_fit_sm_test,Temperature.Ts_fit,Temperature.Patm_final,Temperature.mean_lapse_rate,Temperature.exclusion,Temperature.Titer] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_totalsg,SNRm,cloud_SDm_above);
toc

%%

%=== Cut temperature to surface value ====
Temperature.T_final_test_cut = [Model.Ts(1,:); NaN((cut - 2),Time.i_time); Temperature.T_final_test(cut:end,:)];
Temperature.T_final_test_cut = fillmissing(Temperature.T_final_test_cut,'linear');

gg = cloud_SDm_above.*SNRm;
gg(gg<=0)=nan;
%Appy cloud mask before smoothing
Temperature.T_final_test_cut(isnan(gg)) = NaN;          % Replace mask with NaNs


%Temperature.T_final_testf(isnan(gg)) = NaN;
%Temperature.T_final_testg(isnan(gg)) = NaN;
%[Ez,Et,minSigz,minSigt] = findMinE(Temperature.T_final_testf,Temperature.T_final_testg,0);
%[Temperature.T_final_tests2] = applyFilter(minSigz,minSigt,Temperature.T_final_test_cut);


% k = ones(4,6)./(4*6);     % Kernel
% k = ones(8,6)./(8*6);     % Kernel
% k = ones(8,8)./(8*8);     % Kernel
k = ones(4,6)./(4*6);     % Kernel
%=== apply mask

%==== Smooth temperature
Temperature.T_final_tests = nanconv(Temperature.T_final_test_cut,k,'edge','nanout');

%Temperature.T_final_tests = Temperature.T_final_tests2;
%%%%%Temperature.T_final_tests = Temperature.T_final_test_cut;

%%%Temperature.T_final_tests = Temperature.T_final_test_cut;
%=== apply mask
Temperature.T_finalm = Temperature.T_final_tests ;
Temperature.T_finalm((isnan(gg))) = NaN;
%%
%====== Deconvolution calculation ==========
% % clear o2onDecon o2offDecon
% % load('APDimpulse.mat','impulse_rm','impulse_counts')
% % impulse_counts2 = interp1(impulse_rm*1,impulse_counts,Range.rm);
% % impulse_counts2 = rmmissing(impulse_counts2);
% % impulse_counts2(1) = impulse_counts2(1);
% % impulse_counts2 = impulse_counts2 / trapz(Range.rangeBin,impulse_counts2);
% % 
% % load('AfterPulse.mat')
% % impulse_counts2=pulseON(2:55)/trapz(Range.rangeBin,pulseON(2:end));
% % 
% % for ii=1:Time.i_time
% %     [o2onDecon(:,ii)]=fdeconv(Counts.o2on(:,ii), impulse_counts2)/length(impulse_counts2);
% %     [o2offDecon(:,ii)]=fdeconv(Counts.o2off(:,ii), impulse_counts2)/length(impulse_counts2);
% %     
% % end
% % 
% % o2onDecon = interp2(Time.ts,Range.rm(1:end-length(impulse_counts2)+1),o2onDecon,Time.ts,Range.rm);
% % o2offDecon = interp2(Time.ts,Range.rm(1:end-length(impulse_counts2)+1),o2offDecon,Time.ts,Range.rm);
% % 
% % 
% %  ln_o2_decon = log((o2onDecon(ind_r_lo,:).*o2offDecon(ind_r_hi,:))./(o2onDecon(ind_r_hi,:).*o2offDecon(ind_r_lo,:))); % Natural log of counts
% % alpha_0_raw_decon = ln_o2_decon./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
% % alpha_0_decon = interp2(Time.ts,Range.rm(ind_r_lo),alpha_0_raw_decon,Time.ts,Range.rm);
% % 
% % %%%alpha_0_decon = fillmissing(alpha_0_decon,'nearest');
% % k = ones(8,4)./(8*4);     % Kernel
% % alpha_0_pad_decon = padarray(alpha_0_decon,[8/2,4/2],'replicate');
% % alpha_0_filt_decon = filter2(k,alpha_0_pad_decon,'valid');
% % alpha_0_filt_decon = interp2(Time.ts-Time.t_step/2,Range.rm-Range.rangeBin/2,alpha_0_filt_decon(1:end-1,1:end-1),Time.ts,Range.rm);
% % %%%%%%%%%%%%%%%
% % alpha_0_decon=real(alpha_0_decon);

%%
% ========== Overlap correction calculation ===========
% load('overlapDiffCorr061621.mat','lnOonOoff')
% ln_o2_corr = log((Counts.o2on(ind_r_lo,:).*(Counts.o2off(ind_r_hi,:).*smoothdata(lnOonOoff(ind_r_hi,:),1)))./(Counts.o2on(ind_r_hi,:).*(Counts.o2off(ind_r_lo,:).*smoothdata(lnOonOoff(ind_r_lo,:),1)))); % Natural log of counts
% alpha_0_raw_decon = ln_o2_corr./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
% alpha_0_corr = interp2(Time.ts,Range.rm(ind_r_lo),alpha_0_raw_decon,Time.ts,Range.rm);
% 
% %%%alpha_0_decon = fillmissing(alpha_0_decon,'nearest');
% k = ones(8,4)./(8*4);     % Kernel
% alpha_0_pad_corr = padarray(alpha_0_corr,[8/2,4/2],'replicate');
% alpha_0_filt_corr = filter2(k,alpha_0_pad_corr,'valid');
% alpha_0_filt_corr = interp2(Time.ts-Time.t_step/2,Range.rm-Range.rangeBin/2,alpha_0_filt_corr(1:end-1,1:end-1),Time.ts,Range.rm);
% %%%%%%%%%%%%%%%
% alpha_0_corr=real(alpha_0_corr);

%== Overlap calcuation
for ii=1:size(Sonde.sonde_ind,2)
Sonde.O_on_O_off(:,ii) = Counts.o2on(:,Sonde.sonde_ind(1,ii))./Counts.o2off(:,Sonde.sonde_ind(1,ii)) ./ Sonde.trasmission_sonde{ii}.^2;
end

Model.O_on_O_off = Counts.o2on./Counts.o2off ./ Model.transmission.^2;

%  tic
% %p_point = 481;
% lapse = -11:1:-3; %lapse rate (K) to iterate over
% Stemp = zeros(length(Model.Ts),21);
% Spress = zeros(length(Model.Ts),5);
% for iii=1:length(Model.Ts)
%     Stemp(iii,:) = Model.Ts(iii)-10:1:Model.Ts(iii)+10; %surface temp (K) to iterate over
%     Spress(iii,:) =  Model.Ps(iii)-0.01:0.005:Model.Ps(iii)+0.01; %surface temp (K) to iterate over
% end
% Stemp = fillmissing(Stemp','linear');
% Spress = fillmissing(Spress','linear');
% 
% [~,lowerAlt] = min(abs(1500-Range.rm));
% [~,lowerAlt] = min(abs(1000-Range.rm));
% [~,upperAlt] = min(abs(4000-Range.rm));
% 
% for ii = 1:length(lapse)
%     for jj = 1:size(Stemp,1)
%         for kk = 1:size(Spress,1)
%            TT = Stemp(jj,:) + lapse(ii)/1000 .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
%            PP = Spress(kk,:) .* (Stemp(jj,:)./TT).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
%            alpha = absorption_O2_770_model(TT(1:upperAlt,:),PP(1:upperAlt,:),Spectrum.nu_online(1),Model.WV(1:upperAlt,:)); %[m-1] Funcrtion to calculate theoretical absorption
%            %alpha = abosorption_02_770_PCA(TT(1:upperAlt,:),PP(1:upperAlt,:),Spectrum.nu_online(1),Model.WV(1:upperAlt,:)); %[m-1] Funcrtion to calculate theoretical absorption
%            %[absorption_f,cross_section] = absorption_O2_770_PCA(TT(1:upperAlt,:),TT(1:upperAlt,:),Spectrum.nu_scan_3D_short(1,1,:),Model.WV(1:upperAlt,:));
%            %alpha = absorption_f(:,:,Spectrum.online_index);
%            Transmission = exp(-cumtrapz(Range.rm(1:upperAlt),alpha));
%            dT = diff(log(Transmission.^2),1);
%            measured = log(Range.rm(1:upperAlt).^2 .* Counts.o2on(1:upperAlt,:))-log(Range.rm(1:upperAlt).^2 .* Counts.o2off(1:upperAlt,:));
%            dM = diff(measured,1);
%            
%            meanDiff = sum(abs(dT(1:upperAlt-lowerAlt,:)-dM(1:upperAlt-lowerAlt,:)),1);
%            
%            for iii=1:length(Model.Ts)
%                if ii==1
%                    minmeanDiff(1,iii) = meanDiff(1,iii);
%                    minII(1,iii) = ii;
%                    minJJ(1,iii) = jj;
%                    minKK(1,iii) = kk;
%                elseif meanDiff(1,iii) < minmeanDiff(1,iii)
%                    minmeanDiff(1,iii) = meanDiff(1,iii);
%                    minII(1,iii) = ii;
%                    minJJ(1,iii) = jj;
%                    minKK(1,iii) = kk;
%                end
%            end
%         end
%     end
% end
% 
% %%
% [~,constAlt]=min(abs(2000-Range.rm));
% clear TT PP alpha Transmission 
% TT = zeros(size(Counts.o2on));
% PP = zeros(size(Counts.o2on));
% alpha = zeros(size(Counts.o2on));
% Transmission = zeros(size(Counts.o2on));
% const = zeros(1,size(Counts.o2on,2));
% for iii=1:size(Model.Ts,2)
%     TT(:,iii) = Stemp(minJJ(iii),iii) + lapse(minII(iii))/1000 .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
%     PP(:,iii) = Spress(minKK(iii),iii) .* (Stemp(minJJ(iii),iii)./TT(:,iii)).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
%     alpha(:,iii) = absorption_O2_770_model(TT(:,iii),PP(:,iii),Spectrum.nu_online(1),Model.WV(:,iii)); %[m-1] Funcrtion to calculate theoretical absorption
%     Transmission(:,iii) = exp(-cumtrapz(Range.rm,alpha(:,iii)));
%     const(iii) = log(Range.rm(constAlt).^2 .* Counts.o2on(constAlt,iii))-log(Range.rm(constAlt).^2 .* Counts.o2off(constAlt,iii))-log(Transmission(constAlt,iii).^2);
% end
% 
% overlapCorr = smoothdata(exp(log(Range.rm.^2 .* Counts.o2on)-log(Range.rm.^2 .* Counts.o2off)-log(Transmission.^2)-const),1,'movmean',5);
% 
% toc

%%%%%%%%%end
%%
%===============
%=== Figures ===
%===============

%=Tick spacing and formating
tickHours =4;
tickMin = 0;
[y,m,d]=ymd(span_days(1));
date2=datetime(y,m,d,tickHours,tickMin,0);
date1=datetime(y,m,d,0,0,0);
Format.tickSpacing = datenum(date2)-datenum(date1);
Format.tickAngle = 30;
Format.dateTickFormat ='mm/dd HH:MM';
%%%Format.dateTickFormat ='mm/dd';

%= Plot time for profiles
plot_time = datetime(2021,6,20,14,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
[~,p_point] = min(abs(plot_time-Time.date_ts)); % Find closest value to 338min for comparison to other program
p_point(1:length(Range.rm),1)=p_point;

%= Plot time for profiles with sondes
sonde_index = 1;
p_point = Sonde.sonde_ind(:,sonde_index);
%Sonde.WV_sonde = nan(size(Counts.o2on));

plot_O2(p_point,sonde_index,span_days,Sonde,Model,Counts,Range,Time,Options,Temperature,Format,Alpha,cloud_SDm,HSRL,Data,SNRm,cloud_SDm_above,N_wv,N_wv0,N_wvm,N_wv0m,AbsHumm,AbsHum0m)

%%
function [Ez,Et,minSigz,minSigt] = findMinE(f,g,bg)
    disp('creating filter')
    %====find best filter in range====
    filt_size = logspace(-1,1.5,40); %create filter size in terms of grid points
    Ez = ones(1,length(f(1,:)),length(filt_size));
    for jj = 1:length(filt_size) %loop over different filters
        fprintf('Range %f',jj)
        nz = round(4*filt_size(jj)); %number of grid points 
        z = (-nz:nz)'; %filter grid
        kern = exp(-z.^2/filt_size(jj).^2); %gaussian fitler kernel in range
        if length(kern) > 1 %set gaussian filter for edge cases
            if sum(kern) == 0
                [~,it0] = min(abs(z));
                kern(it0) = 1.0;
            end
        else 
            kern = ones(1);
        end
        kern = kern/sum(sum(kern)); %normalize kernel
        fFilt = nanconv(f-bg,kern,'edge','nanout'); %apply filter
        fFilt(fFilt==0)=.001; %avoid inf in log
        Ez(:,:,jj) = sum(fFilt+bg-g.*log(fFilt+bg),1,'omitnan'); %loss function to optimize
    end
    [~,minEind]=min(Ez,[],3); %find minimum of loss function
    minSigz = ones(1,size(f,2));
    for ii = 1:size(f,2)
        minSigz(ii) = filt_size(minEind(ii)); %set filter width to minimum of loss function
    end

    %====find best filter in time====
    filt_size = logspace(-1.5,5,100);
    filt_size = logspace(-1.5,3,100);
    Et = ones(length(f(:,1)),1,length(filt_size));
    for jj = 1:length(filt_size)
        fprintf('time %f',jj)
        nz = round(4*filt_size(jj)); %number of grid points 
        z = (-nz:nz);%filter grid
        %kern = gaussmf(z,[filt_size(jj),0]);%fitler kernel in range
        kern = exp(-z.^2/filt_size(jj).^2);
        if length(kern) > 1
            if sum(kern) == 0
                [~,it0] = min(abs(z));
                kern(it0) = 1.0;
            end
        else 
            kern = ones(1);
        end
        kern = kern/sum(sum(kern));
        %norm = ones(size(f));
        %norm = conv2(norm,kern,'same');
        %%%%fFilt = conv2(f(:,:)-bg,kern,'same')./norm;
        fFilt = nanconv(f-bg,kern,'edge','nanout');
        fFilt(fFilt==0)=.001;%avoid inf in log
        Et(:,:,jj) = sum(fFilt+bg-g.*log(fFilt+bg),2,'omitnan');
    end
    [~,minEind]=min(Et,[],3);
    minSigt = ones(size(f,1),1);
    for ii = 1:size(f,1)
        minSigt(ii) = filt_size(minEind(ii));
    end       
end

function [counts] = applyFilter(rangeWidth,timeWidth,counts)
    %==Filter in range==
    for iii = 1:size(counts,2) %loop over time
        nz = round(4*rangeWidth(iii)); %number of grid points 
        z = (-nz:nz)';%filter grid
        kern = exp(-z.^2/rangeWidth(iii).^2); %smoothing kernel
        kern = kern/sum(sum(kern)); %normalized smoothing kernel
        counts(:,iii) = nanconv(counts(:,iii),kern,'edge','nanout'); %apply kernel to count data
    end
    %==Filter in time==
    for iii = 1:size(counts,1) %loop over range
        nz = round(4*timeWidth(iii)); %number of grid points 
        z = (-nz:nz);%filter grid
        kern = exp(-z.^2/timeWidth(iii).^2);
        kern = kern/sum(sum(kern));
        counts(:,iii) = nanconv(counts(:,iii),kern,'edge','nanout');
    end
end