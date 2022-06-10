clear all


%===== Varibles
%=Date
date_start = datetime(2021,3,13,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2021,4,3,'TimeZone','UTC');%yyyy,mm,dd

% date_start = datetime(2021,7,12,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,7,13,'TimeZone','UTC');%yyyy,mm,dd
% 
% date_start = datetime(2021,7,20,'TimeZone','UTC');%yyyy,mm,dd
% date_start = datetime(2021,8,29,'TimeZone','UTC');%yyyy,mm,dd
% %date_start = datetime(2021,9,5,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,9,7,'TimeZone','UTC');%yyyy,mm,dd
% 
% date_start = datetime(2021,8,29,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,8,30,'TimeZone','UTC');%yyyy,mm,dd

% date_start = datetime(2021,5,11,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,5,12,'TimeZone','UTC');%yyyy,mm,dd
% 
% date_start = datetime(2021,6,1,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,6,4,'TimeZone','UTC');%yyyy,mm,dd

% date_start = datetime(2021,5,11,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,6,4,'TimeZone','UTC');%yyyy,mm,dd

% date_start = datetime(2021,6,14,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,6,22,'TimeZone','UTC');%yyyy,mm,dd
% 
% date_start = datetime(2021,6,22,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,6,29,'TimeZone','UTC');%yyyy,mm,dd

% % date_start = datetime(2021,5,16,'TimeZone','UTC');%yyyy,mm,dd
% % date_end = datetime(2021,5,17,'TimeZone','UTC');%yyyy,mm,dd
% % 
% % 
% % date_start = datetime(2021,5,11,'TimeZone','UTC');%yyyy,mm,dd
% % date_end = datetime(2021,5,13,'TimeZone','UTC');%yyyy,mm,dd
% % 
% date_start = datetime(2021,3,10,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,3,10,'TimeZone','UTC');%yyyy,mm,dd

% % % date_start = datetime(2021,6,25,'TimeZone','UTC');%yyyy,mm,dd
% % % date_end = datetime(2021,7,2,'TimeZone','UTC');%yyyy,mm,dd

% date_start = datetime(2021,5,18,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,5,23,'TimeZone','UTC');%yyyy,mm,dd
% 
% date_start = datetime(2021,5,21,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2021,5,24,'TimeZone','UTC');%yyyy,mm,dd



span_days = date_start:date_end;

%=Time and range
Options.intTime = 10;  %[min] Integration time
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

[Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadMSUdata(span_days,Options,Constant);
%  HSRL.BSR = ones(size(Counts.o2on));
%  HSRL.BSR(isnan(Counts.o2on))=nan;


% date_begin = datetime(2019,4,17); date_end   = datetime(2019,4,22);
% span_days = date_begin:date_end;        % Days in set [datetime]
% [Range,Time,Counts,Sonde,Model,Spectrum,BSR,Data] = loadSGPdata(span_days,Options,Constant);

% % 
% date_begin = datetime(2020,8,29); 
% date_end   = datetime(2020,9,8);
% % date_begin = datetime(2021,6,16); 
% % %date_begin = datetime(2021,6,19); 
% % date_end   = datetime(2021,6,29);
% % %date_end   = datetime(2021,6,21);
%  span_days = date_begin:date_end;        % Days in set [datetime]
% [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] =loadBoulderdata(span_days,Options,Constant);
% % [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] =loadBoulderdata2(span_days,Options,Constant);



%%
%======================
%= Cloud and SNR mask =
%======================
disp('Calculating masks')
cloud_p_point = 212.5; %hours
%SNR_threshold = 135;
%%SNR_threshold = 7;
SNR_threshold = 4;

SNR_threshold = 36;
SNR_threshold = 31;
SNR_threshold = 36;

%SNR_threshold = .1;
%SNR_threshold = 0;
SD_threshold = 2*10^9;
%SD_threshold = 4*10^10;
%SD_threshold = 10*10^8;
%SD_threshold = 1*10^9;
%%SD_threshold = 5*10^8;
%SD_threshold = 1*10^8;
%SD_threshold = 5*10^20;

SD_threshold = 2;
SD_threshold = 4;
%SD_threshold = 10^5;
BGmult =1;
%BGmult =1;
lowAlt = 0;
%[SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2(Counts.o2on,Counts.o2off,Range.rm,Time.ts,cloud_p_point,SNR_threshold,SD_threshold,Options.oversample,Options.t_avg,Counts,BGmult);
[SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2(cloud_p_point,SNR_threshold,SD_threshold,Options.oversample,Options.t_avg,Counts,BGmult,Time,Range,HSRL.BSR,lowAlt);
%%
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
%==========================
%= Pertabative absorption =
%==========================
disp('Calculating absorption')
for ii=1
    %if ii==2
    %    ln_o2 = log(((Counts.o2on(ind_r_lo,:)-delta(ind_r_lo,:)).*(Counts.o2off(ind_r_hi,:)-delta(ind_r_lo,:)))./((Counts.o2on(ind_r_hi,:)-delta(ind_r_lo,:)).*(Counts.o2off(ind_r_lo,:)+delta(ind_r_lo,:)))); % Natural log of counts
    %else
        ind_r_lo = 1:Range.i_range-Options.oversample;                                            % High range vector
        ind_r_hi = 1+Options.oversample:Range.i_range;                                            % Low range vector
        ln_o2 = log((Counts.o2on(ind_r_lo,:).*Counts.o2off(ind_r_hi,:))./(Counts.o2on(ind_r_hi,:).*Counts.o2off(ind_r_lo,:))); % Natural log of counts
    
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
% % % % %         
%         \

        %end
        
% === Zeroth Order ===

% Zeroth order term

alpha_0_raw = ln_o2./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
alpha_0 = interp2(Time.ts,Range.rm(ind_r_lo),alpha_0_raw,Time.ts,Range.rm);

%%%alpha_0 = fillmissing(alpha_0,'nearest');
k = ones(4,4)./(4*4);     % Kernel
alpha_0_pad = padarray(alpha_0,[4/2,4/2],'replicate');
alpha_0_filt = filter2(k,alpha_0_pad,'valid');
alpha_0_filt = interp2(Time.ts-Time.t_step/2,Range.rm-Range.rangeBin/2,alpha_0_filt(1:end-1,1:end-1),Time.ts,Range.rm);
%%%%%%%%%%%%%%%
alpha_0=real(alpha_0);
%%%%%%%%%

% %Molecular alpha
% %ln_o2_mol = log((Counts.o2on_mol(ind_r_lo,:).*(Counts.o2off_mol(ind_r_hi,:).*BSR(ind_r_hi,:)))./(Counts.o2on_mol(ind_r_hi,:).*(Counts.o2off_mol(ind_r_lo,:).*BSR(ind_r_lo,:)))); % Natural log of counts
% ln_o2_mol = log((Counts.o2on_mol(ind_r_lo,:).*(Counts.o2off_mol_corrected(ind_r_hi,:)))./(Counts.o2on_mol(ind_r_hi,:).*(Counts.o2off_mol_corrected(ind_r_lo,:)))); % Natural log of counts
% alpha_0_raw_mol = ln_o2_mol./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
% alpha_0_mol = interp2(Time.ts,Range.rm(ind_r_lo),alpha_0_raw_mol,Time.ts,Range.rm);
% 
% alpha_0_mol = fillmissing(alpha_0_mol,'nearest');
% k = ones(8,4)./(8*4);     % Kernel
% alpha_0_pad_mol = padarray(alpha_0_mol,[8/2,4/2],'replicate');
% alpha_0_filt_mol = filter2(k,alpha_0_pad_mol,'valid');
% alpha_0_filt_mol = interp2(Time.ts-Time.t_step/2,Range.rm-Range.rangeBin/2,alpha_0_filt_mol(1:end-1,1:end-1),Time.ts,Range.rm);
% %%%%%%%%%%%%%%%
% alpha_0_mol=real(alpha_0_mol);

%%%%%alpha_0 = (alpha_0+alpha_0_mol)/2;


load('TransmittanceData.mat')
T_etalon_on = double(interp1(double(OnlineWavelength)*10^9,OnlineCombinedTransmittance,Spectrum.lambda_scan_3D_short));
T_etalon_off = double(interp1(double(OfflineWavelength)*10^9,OfflineCombinedTransmittance,Spectrum.lambda_scan_3D_short_off));

 altitude = 1.5;%altitude in km
 %[alpha_1, alpha_2] = pertAbsorption(alpha_0, T_etalon_on, T, P, rm,rkm,m_air, nu_online, nu_scan_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,c,kb,altitude);
[alpha_1, alpha_2, Spectrum] = pertAbsorption(alpha_0, T_etalon_on, Model.T, Model.P, Range.rm,Time.ts,Range.rm/1000,Constant.m_air, Spectrum.nu_online, Spectrum.nu_scan_3D_short,Spectrum.nuBin,HSRL.BSR,ind_r_lo,ind_r_hi,Model.WV,Spectrum.online_index,Range.i_range,Time.i_time,Spectrum.i_scan_3D_short,Range.rangeBin,Options.oversample,Options.t_avg,Constant.c,Constant.kb,altitude,Spectrum);

% === Total alpha ===
alpha_total_raw = alpha_0 + alpha_1 + alpha_2;


% Force total alpha to its modeled surface value
[~,cut] = min(abs(Range.rm-1000));             % Index where rm is closest to chosen value
cut=2;
alpha_total_cut = [Model.absorption(1,:); NaN((cut - 2),Time.i_time); alpha_total_raw(cut:end,:)];
alpha_total_cut = fillmissing(alpha_total_cut,'linear');
%alpha_total = alpha_total(2:end,:);     % Remove surface value since rm starts at del_r


% Moving average
% alpha_total_pad = padarray(alpha_total,[oversample/2,t_avg/2],'replicate');
% alpha_total_filt = filter2(k,alpha_total_pad,'valid');
% alpha_total = interp2(ts-t_step/2,rm-rangeBin/2,alpha_total_filt(1:end-1,1:end-1),ts,rm);
% alpha_total = fillmissing(alpha_total,'nearest',1); % Fill in NaNs in dimension 1
% alpha_total = fillmissing(alpha_total,'nearest',2); % Fill in NaNs in dimension 2

%only do zeroth order
%alpha_total = alpha_0;

% ----- smoothing alpha total -----
% for i = 1:i_time
%     alpha_total_pad2(:,i) = padarray(alpha_total(:,i),[4*oversample+1,0],'replicate');
%     alpha_total_pad2(:,i) = smooth(alpha_total_pad2(:,i), 4*oversample+1);
%     alpha_total(:,i) = interp1(rm-rangeBin/2,alpha_total_pad2(4*oversample+1:end-(4*oversample+1)-1,i),rm);
% end

% for i = 1:i_time
%     alpha_total(:,i) = smooth(alpha_total(:,i), 4*oversample+1);
% end

%%

%Smoothing = repmat(gaussmf(linspace(-1,1,8)',[0.5,0]),1,4);
%Smoothing = Smoothing./sum(sum(Smoothing));

%k=Smoothing;

%k = ones(16,8)./(16*8);     % Kernel

k = ones(4,4)./(4*4);     % Kernel

% % % alpha_total_pad = padarray(alpha_total,[6/2,4/2],'replicate');
% % % alpha_total_filt = filter2(k,alpha_total_pad,'valid');
% % % alpha_totals = interp2(Time.ts-Time.t_step/2,Range.rm-Range.rangeBin/2,alpha_total_filt(1:end-1,1:end-1),Time.ts,Range.rm);

alpha_total_cut = alpha_total_cut .* cloud_SDm_above .* SNRm; %.* cloud_SDm_above;
alpha_total_cut(alpha_total_cut == 0) = NaN;          % Replace mask with NaNs

alpha_total_filt = filter2(k,alpha_total_cut,'same');
alpha_totals = alpha_total_filt;


% alpha_total_filt = filter2(k,alpha_total,'same');
% alpha_total = alpha_total_filt;
%T_final_tests = fillmissing(T_final_tests,'nearest',1); % Fill in NaNs in dimension 1
%T_final_tests = fillmissing(T_final_tests,'nearest',2); % Fill in NaNs in dimension 2

%%
%apply SNR mask
alpha_0m = alpha_0 .* SNRm;
alpha_0m(alpha_0m == 0) = NaN;                  % Replace mask with NaNs

alpha_1 = alpha_1 .* SNRm;
alpha_1(alpha_1 == 0) = NaN;                    % Replace mask with NaNs

%alpha_total = real(alpha_total);
alpha_total = real(alpha_totals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_totalm = alpha_total .* cloud_SDm_above .* SNRm; %.* cloud_SDm_above;
alpha_totalm(alpha_totalm == 0) = NaN;          % Replace mask with NaNs

%alpha_total_raw = alpha_total_raw .* SNRm;
%alpha_total_raw(alpha_total_raw == 0) = NaN;    % Replace mask with NaNs

% % % % % for i = 1:i_time
% % % % %     alpha_totalm(:,i) = smooth(alpha_totalm(:,i), 4*oversample+1);
% % % % % end


% % for i = 1:i_time
% %     alpha_totalm(:,i) = smooth(alpha_totalm(:,i), 4*oversample+1);
% % end

%alpha_totalm = alpha_totalm .* SNRm .* cloud_SDm_above;
%alpha_totalm(alpha_totalm <= 0) = NaN;          % Replace mask with NaNs

%%

%delta = findDelta(Model.absorption(:,:),alpha_total(:,:),Counts.o2on(:,:),Counts.o2off(:,:),Range.rm,Time.ts);
delta = findDelta(Model.absorption(:,:),alpha_total_raw(:,:),Counts.o2on(:,:),Counts.o2off(:,:),Range.rm,Time.ts);
delta=-delta;

delta = [delta; zeros(length(Range.rm)-size(delta,1),length(Time.ts))];
delta = delta/4;


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
end

%%
disp('Temp retrieval function')
%WV = 0; % set water vapor concentration to zero
tic
[T_final_test,L_fit_sm_test,Ts_fit,Patm_final,mean_lapse_rate,exclusion] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,alpha_totalm,SNRm,cloud_SDm_above);
toc



% % % alpha_total_pad = padarray(alpha_total,[oversample/2,t_avg/2],'replicate');
% % % alpha_total_filt = filter2(k,alpha_total_pad,'valid');
% % % alpha_total = interp2(ts-t_step/2,rm-rangeBin/2,alpha_total_filt(1:end-1,1:end-1),ts,rm);
% % % alpha_total = fillmissing(alpha_total,'nearest',1); % Fill in NaNs in dimension 1
% % % alpha_total = fillmissing(alpha_total,'nearest',2); % Fill in NaNs in dimension 2
% 
% T_final_test_pad = padarray(T_final_test,[8/2,30/2],'replicate');
% T_final_test_filt = filter2(k,T_final_test_pad,'valid');
% T_final_test_filt = interp2(ts-t_step/2,rm-rangeBin/2,T_final_test_filt(1:end-1,1:end-1),ts,rm);
% T_finalm = T_final_test_filt .* SNRm .* cloud_SDm_above;
% T_finalm(T_finalm <= 0) = NaN;

T_finalm = T_final_test .* SNRm .* cloud_SDm_above;
T_finalm(T_finalm <= 0) = NaN;


%k = ones(14,10)./(14*10);     % Kernel

%Smoothing = repmat(gaussmf(linspace(-1,1,8)',[0.5,0]),1,4);
%Smoothing = repmat(gaussmf(linspace(-1,1,14)',[0.5,0]),1,10);
%Smoothing = Smoothing./sum(sum(Smoothing));

%k=Smoothing;

T_final_test_cut = [Model.Ts(1,:); NaN((cut - 2),Time.i_time); T_final_test(cut:end,:)];
T_final_test_cut = fillmissing(T_final_test_cut,'linear');

%k = ones(8,4)./(8*4);     % Kernel
k = ones(4,4)./(4*4);     % Kernel
% 
% %T_final_test_pad = padarray(T_final_test,[14/2,10/2],'replicate');
% T_final_test_pad = padarray(T_final_test,[8/2,8/2],'replicate');
% T_final_test_filt = filter2(k,T_final_test_pad,'valid');
% T_final_tests = interp2(Time.ts-Time.t_step/2,Range.rm-Range.rangeBin/2,T_final_test_filt(1:end-1,1:end-1),Time.ts,Range.rm);
% %%%T_final_tests = fillmissing(T_final_tests,'nearest',1); % Fill in NaNs in dimension 1
% %%%T_final_tests = fillmissing(T_final_tests,'nearest',2); % Fill in NaNs in dimension 2


T_final_tests = filter2(k,T_final_test_cut,'same');


T_finalm = T_final_tests .* SNRm .* cloud_SDm_above ;

%T_finalm = T_final_test .* SNRm .* cloud_SDm_above ;
T_finalm(T_finalm <= 0) = NaN;




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
%%
ts = Time.ts;
thr = Time.ts/60/60;
rm =Range.rm;
rkm = rm/1000;

reportfile = 'C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\5_7_21\';

%Tick spacing
tickHours = 8;
tickMin = 0;
[y,m,d]=ymd(span_days(1));
date2=datetime(y,m,d,tickHours,tickMin,0);
date1=datetime(y,m,d,0,0,0);
tickSpacing = datenum(date2)-datenum(date1);
tickAngle = 30;
dateTickFormat ='mm/dd HH:MM';

%===============
%=== Figures ===
%===============

plot_time=5;
[~,p_point] = min(abs(plot_time-Time.ts/60/60)); % Find closest value to 338min for comparison to other program
p_point(1:length(Range.rm),1)=p_point;

plot_time = datetime(2021,6,23,8,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2021,8,6,8,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2020,9,4,8,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2021,8,11,20,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2021,8,30,1,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
%plot_time = datetime(2021,8,4,2,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
[~,p_point] = min(abs(plot_time-Time.date_ts)); % Find closest value to 338min for comparison to other program
p_point(1:length(Range.rm),1)=p_point;

sonde_index = 1;
%p_point = Sonde.sonde_ind(:,sonde_index);
%Sonde.sonde_ind = [];

%==Fitted lapse rate and surface temperature
figure(728590)
subplot(2,1,1)
plot(ts/60/60,permute(L_fit_sm_test(1,:,end)*1000,[3 2 1]))
hold on
plot(ts([p_point(1) p_point(end)])/60/60,[-10 -4],'k','linewidth',2)
hold off
yline(-6.5);
grid on
title('Fitted lapse rate')
title(sprintf('Fitted lapse rate\n %s',span_days(1)))
xlabel('Time (UTC hr)')
ylabel('Lapse rate (K/km)')
subplot(2,1,2)
plot(ts/60/60,Model.T(1,:) - permute(Ts_fit(1,:,:),[3 2 1]))
hold on
plot(ts([p_point(1) p_point(end)])/60/60,[-20 20])
hold off
grid on
title('Surface temp - fitted surface temp')
xlabel('Time (UTC hr)')
ylabel('\deltaT (K)')

%==Background Counts
figure(9985)
semilogy(datenum(Time.date_ts),Counts.bg_o2off,'-',datenum(Time.date_ts),Counts.bg_o2on,datenum(Time.date_ts),Counts.bg_o2off_mol,datenum(Time.date_ts),Counts.bg_o2on_mol,datenum(Time.date_ts([p_point(1) p_point(end)])),[1 100])
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
xline(738397.94)
xline(738398.17)
datetick('x',dateTickFormat,'keeplimits','keepticks')
title('Background')
xlabel('Time UTC hr')
ylabel('photons')
legend('Offline','Online','Offline molecular','Online molecular')
grid on

%=Diff Counts and sonde difference
figure(99496)
cla reset
line(diff(Counts.o2off(:,p_point(1)))/Range.rangeBin,rkm(1:end-1),'Color','r')
line(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin./Counts.o2off(1:end-1,p_point(1))*100,rkm(1:end-1),'Color','b')
%line(diag(dBSRSG(:,p_point)./BSR(:,p_point)),rkm,'Color','b')
%line(diag(dBSRlo(:,p_point)),rkm(1:end-1))
legend('dOFF/dr')
xlim([-2 2])
xlabel('dOFF/dr')
ylabel('Range (km)')
ax1 = gca; % current axes
ax1.XColor = 'r';
ax1.YColor = 'r';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
%line(diag(alpha_0(:,p_point)-Sonde.absorption_sonde{sonde_index}),rkm,'Parent',ax2,'Color','k')
line(diag(alpha_0_filt(:,p_point)-Sonde.absorption_sonde{sonde_index}),rkm,'Parent',ax2,'Color','g')
%line((diag(alpha_0_filt(1:end-1,p_point)))./(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin./Counts.o2off(1:end-1,p_point(1))*100)-Sonde.absorption_sonde{sonde_index}(1:end-1),rkm(1:end-1),'Parent',ax2,'Color','k')
%line((diag(alpha_0_filt(1:end-1,p_point)))./(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin)-Sonde.absorption_sonde{sonde_index}(1:end-1),rkm(1:end-1),'Parent',ax2,'Color','k')
grid on
xlim([-.5 .5]*10^-3)
xlabel('\alpha_0 - \alpha_{sonde}')
ylabel('Range (km)')
legend('\alpha_0 - \alpha_{sonde}')

%=Water Vapor profile
figure(99499)
subplot(2,1,1)
plot(diag(Model.WV(:,p_point)),rm)
xlabel('WV (molec/m^3')
ylabel('Range (m)')
legend('VW')
xlabel('Absolute Humidity g/m^3')

%=Temperature profile
figure(884)
%%%plot(diag(Model.T(:,p_point)),rkm)
%plot(T(:,p_point)+2,rkm)
%plot(T(:,p_point)-2,rkm)
%plot(Ts_fit(1,p_point,end),0,'+')
%plot(diag(T_final_test(:,p_point)),rkm)
plot(diag(T_finalm(:,p_point)),rkm,'linewidth',2)
exclusion(exclusion==0)=NaN;
%plot(diag(T_final_test(:,p_point)).*diag(exclusion(:,p_point,end)),rkm,'*')
%plot(Sonde.T_sonde(:,sonde_index),rm/1000,'.-')
%plot(Model.Ts(p_point(1)),rkm(1),'+')
hold on
plot(Model.T(:,p_point(1)),rkm,'-','LineWidth',2)
plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end)),rkm,'--')
hold off
grid on
xlabel('Temperature (K)')
ylabel('Range (km)')
ylim([0 5])
title(sprintf(['Temp guess and retrieved from pertabative absorption\n' datestr(Time.date_ts(p_point(1)))]))
legend('Retrieved Temperature','Model temperature','Location','northeast')
%%%saveas(gcf,[reportfile sprintf('temperature%.0f0%.0f.png',floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60)])

%=Tfit-T
figure(885)
plot(diag(L_fit_sm_test(:,p_point,end) .* rm + Ts_fit(1,p_point,end))-diag(T_finalm(:,p_point)),rkm)
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-1 -1],[0 6],'Color','k','LineStyle','-.')
line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color','k')
line([2 2],[0 6],'Color','k')
xlabel('\Delta T (T_{fit} - T_{retrieved}) (K)')
ylabel('Range (km)')
legend('\Delta T','\Delta T = 0','\Delta T = \pm 1','\Delta T = \pm 2','Location','northwest')
ylim([0 6])
xlim([-10 10])
title(sprintf(['Temp fit and DIAL difference\n' datestr(Time.date_ts(p_point(1)))]))
hold off

%=Tsonde-T
figure(888)
plot(Sonde.T_sonde(:,sonde_index) - diag(T_finalm(:,p_point)),rkm)
plot(diag(Model.T(:,p_point)) - diag(T_finalm(:,p_point)),rkm)
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-1 -1],[0 6],'Color','k','LineStyle','-.')
line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color','k')
line([2 2],[0 6],'Color','k')
xlabel('\Delta T (T_{sonde} - T_{retrieved}) (K)')
ylabel('Range (km)')
ylim([0 5])
xlim([-10 10])
title(sprintf(['Temperature difference\n' datestr(Time.date_ts(p_point(1)))]))
%%%saveas(gcf,[reportfile sprintf('tempDiff%.0f0%.0f.png',floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60)])

%=2D temperature
figure(886)
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot)); % Max range Index
r_min_plot = .5; % Max range to plot [km]
[~,ind_km_min] = min(abs(rkm-r_min_plot)); % Max range Index
imAlpha=ones(size(T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
imAlpha(isnan(T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
%imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
colorLimits = [min(T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1
colorLimits =[-25 50]+273;
%colorLimits =[-5 30]+273;
imagesc(datenum(Time.date_ts),rkm(1:ind_km_max),T_finalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits) %Plot
hold on
%sonde lines
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[rkm(1) rkm(end)],'--b','LineWidth',2)
%%%%%%%%%%%%%%%%%%plot(datenum(Time.date_ts([833 833])),[rkm(1) rkm(end)],'--k','LineWidth',2)
% % yline(.6)
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[rkm(1) rkm(end)],'--k')
% % end

xline(738397.94)
xline(738398.17)

xline(738400.65)
xline(738401.18)
colors = colormap; % Set colors to colormap
colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
colormap(colors); % Make colormap new colors
set(gca, 'YDir','normal') % Set y axis increasing
set(gca,'color',[1 1 1]);% Color background white
colorbar % Add colorbar
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
%title(sprintf('Temperature retrieval\n%s(UTC)',span_days(1))) %Add title
title(sprintf('Temperature retrieval')) %Add title
xlabel('Time (UTC hours)')
ylabel('Range (km)')
title(colorbar,'Temperature (K)') % Add title to colorbar
hold off

%=Absorption profile
figure(4)
plot(diag(alpha_0(:,p_point)),rkm)
hold on
plot(diag(alpha_totalm(:,p_point)),rkm,'LineWidth',2)
%plot(diag(alpha_total_raw(:,p_point)),rkm)
%%%plot(Sonde.absorption_sonde{sonde_index},rkm,'.-')
%plot(diag(Model.absorption(:,p_point)),rkm)
%plot(diag(alpha_d(:,p_point)),rkm,'--','LineWidth',2)
%plot(diag(alpha_0_filt(:,p_point)),rkm,'--','LineWidth',2)
%plot(diag(alpha_0_filt_mol(:,p_point)),rkm,'--','LineWidth',2)
%plot(alpha_0_filt_decon(:,p_point(1)),rkm,'.-')
%%%plot(alpha_0_filt_corr(:,p_point(1)),rkm,'k.-','LineWidth',2)
plot(Model.absorption(:,p_point(1)),rkm,'-','LineWidth',2)
%plot(alpha_d(:,p_point(1)),rkm(1:end-1),'.-')
%plot(alpha_dModel(:,p_point(1)),rkm(1:end-1),'.-')
hold off
grid on
title(sprintf(['Calculated O_2 absorption coefficient\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')
xlim([-2e-4 6e-4])
ylim([0 5])
%legend('Measured zeroth order absorption','Pertabative absorption','Alpha 0+1+2','Radiosonde model','\alpha model','\alpha_0+\Delta','Smooth \alpha_0')
legend('Measured zeroth order absorption','Purtabative absorption','Model absorption')
%%%saveas(gcf,[reportfile sprintf('absorption%.0f0%.0f.png',floor(thr(p_point(1))),(thr(p_point(1))-floor(thr(p_point(1))))*60)])

figure(44)
plot((Model.absorption(:,p_point(1))-alpha_totalm(:,p_point(1)))./Model.absorption(:,p_point(1))*100,rkm,'-')
hold on
plot((Model.absorption(:,p_point(1))-alpha_0(:,p_point(1)))./Model.absorption(:,p_point(1))*100,rkm,'-')
%%plot((Model.absorption(1:end-1,p_point(1))-alpha_d(:,p_point(1)))./Model.absorption(1:end-1,p_point(1))*100,rkm(1:end-1),'-')
hold off
legend('model total diff','model 0 diff','d')
xlabel('percent diff from model')
grid on
xlim([-200 200])


%=BSR
figure(444)
plot(diag(HSRL.BSR(:,p_point)),rkm,'linewidth',2)
hold off
grid on
title(sprintf(['BSR\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlim([0 10])
ylim([0 5])


%=Counts Profile
figure(7)
subplot(2,2,1)
plot(diag(Counts.o2on(:,p_point)),rm)
hold on
plot(diag(Counts.o2off(:,p_point)),rm)
plot(diag(Counts.o2on_mol(:,p_point)),rm)
plot(diag(Counts.o2off_mol(:,p_point)),rm)
plot(Counts.o2on_bgsub(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2off_bgsub(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2on_bgsub_mol(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2off_bgsub_mol(:,p_point(1)),Range.rm_raw_o2,'--')
xlim([0 500])
%view([90 -90])
title('Averaged photon returns')
legend('Online','Offline','Online molecular','Offline molecular','raw online','raw offline','raw online molecular','raw offline molecular')
hold off
grid on
ylabel('Range [m]')
xlabel('Photons')
subplot(2,2,2)
plot(log(Counts.o2on_bgsub(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
hold on
plot(log(Counts.o2off_bgsub(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
plot(log(Counts.o2on_bgsub_mol(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
plot(log(Counts.o2off_bgsub_mol(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('log(Photons*range^2)')
ylabel('Range [m]')
grid on
hold off
subplot(2,2,3)
semilogx((Counts.o2on_bgsub(:,p_point(1))),Range.rm_raw_o2/1000);
hold on
semilogx((Counts.o2off_bgsub(:,p_point(1))),Range.rm_raw_o2/1000);
semilogx((Counts.o2on_bgsub_mol(:,p_point(1))),Range.rm_raw_o2/1000);
semilogx(Counts.o2off_bgsub_mol(:,p_point(1)),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('Photons')
ylabel('Range [m]')
grid on
hold off
subplot(2,2,4)
semilogx(Counts.o2on_bgsub(:,p_point(1))+Counts.bg_o2on(:,p_point(1)),Range.rm_raw_o2/1000);
hold on
semilogx(Counts.o2off_bgsub(:,p_point(1))+Counts.bg_o2off(:,p_point(1)),Range.rm_raw_o2/1000);
semilogx(Counts.o2on_bgsub_mol(:,p_point(1))+Counts.bg_o2on_mol(:,p_point(1)),Range.rm_raw_o2/1000);
semilogx(Counts.o2off_bgsub_mol(:,p_point(1))+Counts.bg_o2off_mol(:,p_point(1)),Range.rm_raw_o2/1000);

% semilogx(Counts.o2on_raw(:,p_point(1)*10),Range.rm_raw_o2/1000,'--');
% semilogx(Counts.o2off_raw(:,p_point(1)*10),Range.rm_raw_o2/1000,'--');
legend('on','off','on mol','off mol')
xlabel('Photons')
ylabel('Range [m]')
grid on
hold off

%=Housekeeping temperature
figure(738732)
subplot(4,1,1)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.InsideCell.Temperature)
legend('Inside Cell')
title('Thermocouples')
ylabel('^o C')
subplot(4,1,2)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.OutsideCell.Temperature  )
ylabel('^o C')
legend('Outside Cell')
subplot(4,1,3)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.TSOA.Temperature  )
ylabel('^o C')
legend('TSOA')
subplot(4,1,4)
plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.RoomTemp.Temperature  )
ylabel('^o C')
legend('Room Temp')
xlabel('Time UTC hr')

%=Housekeeping temperature
% figure(738732)
% subplot(8,1,1)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.InsideCell.Temperature)
% legend('HSRLOven')
% title('Thermocouples')
% ylabel('^o C')
% subplot(8,1,2)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.HVACReturn.Temperature  )
% ylabel('^o C')
% legend('HVAC return')
% subplot(8,1,3)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.HVACSource.Temperature  )
% ylabel('^o C')
% legend('HVACSource')
% subplot(8,1,4)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.Window.Temperature  )
% ylabel('^o C')
% legend('window')
% subplot(8,1,5)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.WindowHeaterFan.Temperature  )
% ylabel('^o C')
% legend('windowHeaterfan')
% subplot(8,1,6)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.VWEtalonHeatSink.Temperature  )
% ylabel('^o C')
% legend('WVetalon heat sink')
% subplot(8,1,7)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.HSRLEtalonHeatSink.Temperature  )
% ylabel('^o C')
% legend('HSRL etalon heat sink')
% subplot(8,1,8)
% plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.RoomTemp    .Temperature  )
% ylabel('^o C')
% legend('optical table')
% xlabel('Time UTC hr')

%=Laser locking wavelenths
figure(738734)
subplot(2,1,1)
plot(Data.Laser.O2Online.TimeStamp,Data.Laser.O2Online.WaveDiff)
hold on
plot(Data.Laser.O2Online.TimeStamp,Data.Laser.O2Offline.WaveDiff)
hold off
legend('Online','Offline')
title('wavelength Difference')
ylim([-5*10^-4 5*10^-4])
ylabel('(nm)')
xlabel('Time hr')
subplot(2,1,2)
plot(Data.Laser.O2Online.TimeStamp  ,Data.Laser.O2Online.WaveDiff-Data.Laser.O2Offline.WaveDiff  )
title('Online wave diff-offline wave diff')
ylim([-5*10^-4 5*10^-4])
grid on
ylabel('(nm)')
xlabel('Time hr')

%=Laser temperatures
figure(738735)
subplot(4,1,1)
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Online.TemperatureActual)
hold on
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Online.TemperatureDesired)
legend('Actual','Desired')
hold off
grid on
title('Online Laser Temperature')
ylabel('(C)')
subplot(4,1,2)
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Offline.TemperatureActual)
hold on
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Offline.TemperatureDesired)
legend('Actual','Desired')
hold off
grid on
title('Offline Laser Temperature')
ylabel('(C)')
subplot(4,1,3)
% plot(Data.Laser.TWSOA.TimeStamp,Data.Laser.TWSOA.TemperatureActual)
% hold on
% grid on
% plot(Data.Laser.TWSOA.TimeStamp,Data.Laser.TWSOA.TemperatureDesired)
% legend('Actual','Desired')
% hold off
% title('TWSOA Laser Temperature')
% ylabel('(C)')
subplot(4,1,4)
plot(Data.Etalon.O2Etalon.TimeStamp,Data.Etalon.O2Etalon.TemperatureActual)
title('Etalon Temperature')
grid on
ylabel('(C)')
xlabel('Time hr')

%=Laser seed power
figure(778298)
plot(Data.Laser.O2Online.TimeStamp,10.^(Data.Laser.O2Online.SeedPower/10))
hold on
plot(Data.Laser.O2Offline.TimeStamp,10.^(Data.Laser.O2Offline.SeedPower/10))
hold off
%ylim([-5.5 -4.8])
title('Laser Seed power at wavemeter')
legend('Online','Offline')

%=Weather Station
figure(83)
subplot(2,1,1)
plot(datenum(Time.date_ts),Model.Ts)
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
title('Weather Station Surface Temperature')
xlabel('Time')
ylabel('K')
subplot(2,1,2)
plot(datenum(Time.date_ts),Model.Ps)
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
title('Weather Station Suface Pressure')
xlabel('Time')
ylabel('Pressure ATM')

%=UPS
if ishandle(648)
    figure(648)
    clf(648)
end
figure(648)
line(Data.UPS.all.TimeStamp,Data.UPS.all.BatteryTimeLeft,'Color','r')
line(Data.UPS.all.TimeStamp,Data.UPS.all.HoursOnBattery,'Color','g')
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('Battery remaining (hrs)')
legend('Battery remaining','Hours on battery','Location','Southwest')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(Data.UPS.all.TimeStamp,Data.UPS.all.UPSTemperature,'Parent',ax2,'Color','k')
ylabel('Battery temperature (C)')

%=2d delta
figure(4489)
imagesc(thr,rm/1000,delta)
hold on
xline(thr(p_point(1)),'K');
hold off
set(gca, 'YDir','normal')
caxis([-1000 1000]/2)
colorbar
%colormap('hot')
title('delta')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')

%=2D counts
figure(589)
subplot(2,2,1)
imagesc(thr,rm/1000,Counts.o2on_noise.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
colorbar
colormap('hot')
title('o2 on unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(thr,rm/1000,Counts.o2on_noise)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 500])
colorbar
colormap('hot')
title('o2 on unaveraged ')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,3)
imagesc(thr,rm,Counts.o2off_noise.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
colorbar
title('o2off unaveraged  * r^{ 2}')
subplot(2,2,4)
imagesc(thr,rm,Counts.o2off_noise)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 500])
colorbar
title('o2off unaveraged')

%=2D molecular counts
figure(590)
subplot(2,2,1)
imagesc(thr,rm/1000,Counts.o2on_noise_mol.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 1e9])
colorbar
colormap('hot')
title('o2 on mol unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(thr,rm/1000,Counts.o2on_noise_mol)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 500])
colorbar
colormap('hot')
title('o2 on mol unaveraged')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,3)
imagesc(thr,rm,Counts.o2off_noise_mol.*rm.^2)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 5e7])
colorbar
colormap('hot')
title('o2 off mol unaveraged * r^{ 2}')
subplot(2,2,4)
imagesc(thr,rm,Counts.o2off_noise_mol)
hold on
xline(thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
caxis([0 75])
colorbar
colormap('hot')
title('o2 off mol unaveraged')

%=2D absorption
figure(8886)
subplot(3,1,1)
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(rkm-r_max_plot));
imAlpha=ones(size(alpha_0m(1:ind_km_max,:)));
imAlpha(isnan(alpha_0m(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_0m(1:ind_km_max,:),[],'all'), max(alpha_0m(1:ind_km_max,:),[],'all')];
colorLimits = [-3e-3, 3e-3];
imagesc(datenum(Time.date_ts),rkm(1:ind_km_max),alpha_0m(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[rkm(1) rkm(end)],'--b')
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha 0 %s to %s (UTC)',Time.date_ts(1),Time.date_ts(end)))
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xlim([datenum(Time.date_ts(1)) datenum(Time.date_ts(end))]) 
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
hold off
subplot(3,1,2)
alpha_total_rawm = alpha_total_raw .* SNRm .* cloud_SDm_above;
imAlpha=ones(size(alpha_total_rawm(1:ind_km_max,:)));
imAlpha(isnan(alpha_total_rawm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_total_rawm(1:ind_km_max,:),[],'all'), max(alpha_total_rawm(1:ind_km_max,:),[],'all')];
colorLimits = [min(alpha_total_rawm(1:ind_km_max,:),[],'all')/5, max(alpha_total_rawm(1:ind_km_max,:),[],'all')/5];
imagesc(datenum(Time.date_ts),rkm(1:ind_km_max),alpha_total_rawm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[rkm(1) rkm(end)],'--b')
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha raw'))
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
hold off
subplot(3,1,3)
imAlpha=ones(size(alpha_totalm(1:ind_km_max,:)));
imAlpha(isnan(alpha_totalm(1:ind_km_max,:)))=0;%Set AlphaData to not plot NaNs
imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1;
colorLimits = [min(alpha_totalm(1:ind_km_max,:),[],'all'), max(alpha_totalm(1:ind_km_max,:),[],'all')];
colorLimits = [-2e-4, 10e-4];
imagesc(datenum(Time.date_ts),rkm(1:ind_km_max),alpha_totalm(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[rkm(1) rkm(end)],'--b')
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
set(gca, 'YDir','normal')
set(gca,'color',[1 1 1]);%color background white
colorbar
title(sprintf('alpha total'))
xlabel('Time UTC')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
hold off

%=2D BSR
figure(7473)
BSRm = HSRL.BSR.*cloud_SDm_above.*SNRm;
imagesc(datenum(Time.date_ts),rkm,BSRm)
hold on
plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[rkm(1) rkm(end)],'--b')
hold off
colorbar
colormap(flipud(hot))
caxis([0 10])
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
datetick('x',dateTickFormat,'keeplimits','keepticks')
ylabel('Range (m)')
set(gca, 'YDir','normal')
title('BSR')

%=T sonde and DIAL profiles
figure(685438)
plot(240:300,240:300,'--k')
hold on
plot(242:302,240:300,'k')
plot(238:298,240:300,'k')
for ii = 1:size(Sonde.sonde_ind,2)
    plot(diag(T_finalm(:,Sonde.sonde_ind(:,ii))),Sonde.T_sonde(:,ii))
end
hold off
xlabel('T DIAL')
ylabel('T Sonde')

%=T Sonde and DIAL scatter plot
figure(67755)
[~,sondeCut]=min(abs(rm-500));
tempProb = zeros(330-240+1);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:length(240:330)
        for j = 1:length(240:330)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(240+i-1)+1)&(diag(T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(240+i-1)-0)& (Sonde.T_sonde(sondeCut:end,ii)<(240+j-1)+1)&(Sonde.T_sonde(sondeCut:end,ii)>=(240+j-1)-0));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(240:330,240:330,tempProb,'AlphaData',imAlpha)
colormap(flipud(hot))%colormap hot
set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
plot(239:302,239:302,'k','LineWidth',2)
plot(240:303,237:300,'k--','LineWidth',2)
plot(237:300,240:303,'k--','LineWidth',2)
hold off
xlim([239 302])
ylim([239 302])
set(gca, 'YDir','normal')
xlabel('Temp retrieval (K)')
ylabel('Sonde (K)')

%=T Sonde-DIAL histogram
figure(67756)
hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
for ii = 1:size(Sonde.sonde_ind,2)
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Sonde.T_sonde(jj,ii)-T_finalm(jj,Sonde.sonde_ind(jj,ii));
        int=int+1;
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('Histogram T_{sonde}-T_{DIAL}\n Mean %f, std %f',meanHist,stdHist))
end

%=FFT of Counts
figure(7726)
subplot(2,2,1)
Fs = 1/Range.rangeBin;
T = Range.rangeBin;
L = length(Range.rm);
f = Fs*(0:(L/2))/L;
Y = fft(Counts.o2on(:,p_point(1)).*rm.^2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold on
Y = fft(Counts.o2off(:,p_point(1)).*rm.^2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
Y = fft(Counts.o2on_mol(:,p_point(1)).*rm.^2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
Y = fft(Counts.o2off_mol(:,p_point(1)).*rm.^2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold off
title('FFT of counts in range at plot point')
legend('on','off','on mol','off mol')
xlabel('f (1/m)')
ylabel('|P1(f)|')
%figure(722897)
subplot(2,2,2)
Y = fft(fillmissing(alpha_total(:,p_point(1)),'nearest'));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold on
Y = fft(fillmissing(alpha_total_raw(:,p_point(1)),'nearest'));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
Y = fft(fillmissing(Sonde.absorption_sonde{sonde_index}(:,:),'nearest'));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold off
legend('Alpha total','alpha total raw','alpha sonde')
title('FFT of absorption in range dimension at plot point')
xlabel('Freqency (1/m)')
%figure(722898)
subplot(2,2,3)
Y = fft(fillmissing(T_final_tests(:,p_point(1)),'nearest'));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold on
Y = fft(fillmissing(T_final_test(:,p_point(1)),'nearest'));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
Y = fft(fillmissing(Sonde.T_sonde(:,sonde_index),'nearest'));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold off
legend('T total','T total raw','T sonde')
title('FFT of Temperature in range dimension at plot point')
xlabel('Freqency (1/m)')
%figure(722899)
subplot(2,2,4)
Y = fft([ones(1,8) zeros(1,Range.i_range-8)]);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold on
Y = fft([ones(1,20) zeros(1,Range.i_range-20)]);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
semilogy(f,P1)
hold off
legend('8 length','10 length')

%=N bins over time
figure(66572)
bar(datenum(Time.date_ts),Counts.NBins)
ylabel('# of bins summed')
xticks(datenum(Time.date_ts(1)): tickSpacing :datenum(Time.date_ts(end)))
xtickangle(tickAngle)
xline(738397.94)
xline(738398.17)

xline(738400.65)
xline(738401.18)
datetick('x',dateTickFormat,'keeplimits','keepticks')

%=Counts and corrections
figure(6)
plot(Range.rm,diag(Counts.o2on(:,p_point)))
hold on
plot(Range.rm,diag(Counts.o2off(:,p_point)))
plot(Range.rm,diag(Counts.o2on_mol(:,p_point)))
plot(Range.rm,diag(Counts.o2off_mol(:,p_point)))
%plot(Range.rm,diag(delta(:,p_point)),'.-')
%plot(Range.rm,diag(deltaModel(:,p_point)),'.-')
%plot(Range.rm,diag(Counts.o2on(:,p_point)-delta(:,p_point)),'--')
%plot(Range.rm,diag(Counts.o2off(:,p_point)-delta(:,p_point)),'--')
ylim([-200 1000])
%load('AfterPulse.mat')
%plot(rm,pulseON(4:160+3))
%plot(rm,[pulseON(4:160); zeros(3,1)])
% plot(Range.rm,diag(Counts.o2on(:,p_point)-pulseON(4:133+3)),'.-')
% plot(Range.rm,diag(Counts.o2off(:,p_point)-pulseON(4:133+3)),'.-')
%%plot(Range.rm,diag(Counts.o2on(:,p_point)-[pulseON(4:160); zeros(3,1)]),'.-')
%%plot(Range.rm,diag(Counts.o2off(:,p_point)-[pulseON(4:160); zeros(3,1)]),'.-')
% plot(Range.rm,diag(Counts.o2on_mol(:,p_point)))
% plot(Range.rm,diag(Counts.o2off_mol(:,p_point)))
%plot(Range.rm(1:end-1),-diff(Counts.o2off(:,p_point(1))))
%%plot(Range.rm(1:end-1),-diff(Counts.o2off(:,p_point(1)))-[pulseON(4:159); zeros(3,1)])

plot(Range.rm,Model.N_on(:,p_point(1)))
plot(Range.rm,Model.N_off(:,p_point(1)))
plot(Range.rm,Model.N_on_pulse(:,p_point(1)),'--')
plot(Range.rm,Model.N_off_pulse(:,p_point(1)),'--')
grid on
hold off
%legend('On','off','\Delta','On -\Delta','Off -\Delta','Measured afterpulse')
legend('On','off','on model','off model','on model pulse','off model pulse')
title(sprintf(['Counts\n' datestr(Time.date_ts(p_point(1)))]))

figure(666)
plot(Range.rm,Model.OverlapOn(:,p_point(1)));
hold on
plot(Range.rm,Model.OverlapOff(:,p_point(1)));
plot(Range.rm,Model.OverlapOn_pulse(:,p_point(1)),'--');
plot(Range.rm,Model.OverlapOff_pulse(:,p_point(1)),'--');
hold off
title('Model counts overlap')
legend('on','off','on pulse','off pulse')

%=Wavelength overlap correction model
% figure(555)
% plot(rkm,Model.O_on_O_off(:,p_point(1)-10:p_point(1)+10))
% grid on
% title('O_{on}(r)/O_{off}(r) correction')
% xlabel('Range km')

%=Transmission
figure(90)
lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(Model.transmission.^2);
plot(Range.rm,lnOonOoff(:,p_point(1)))
hold on
for time = 1:length(ts)
        % Set fit exclusion zones
        exclusion(:,time) = rm<1500 | rm>2500 | SNRm(:,time)==0 | cloud_SDm_above(:,time) ~= 1;
        logicalExc(:,time) = ~logical(exclusion(:,time)); 
        V = [rm(logicalExc(:,time)) ones(length(rm(logicalExc(:,time))),1)];
        %calculate polynomial coefficients     
        p = V\alpha_0(logicalExc(:,time),time);
        alpha_fit_slope(:,time) = p(1);
        alpha_fit_int(:,time) = p(2);
        alpha_fit(:,time) = alpha_fit_int(:,time) + rm.*alpha_fit_slope(:,time);
end
transmission = exp(-cumtrapz(Range.rm,alpha_fit));
lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(transmission.^2);
%lnOonOoff = lnOonOoff - mean(lnOonOoff(40:67,:),1);
lnOonOoff = lnOonOoff - mean(lnOonOoff(14:27,:),1);
plot(Range.rm,lnOonOoff(:,p_point(1)),'--')
grid on
legend('model lnOonOoff','Fit lnOonOoff')
hold off

%=log on and off
figure(91)
lnOonOoff=exp(lnOonOoff);
plot(Range.rm,lnOonOoff(:,p_point),'--')
%plot(Range.rm,lnOonOoff(:,p_point)./Corr(:,p_point))
title('lnOonOoff')
hold off
           
% figure
% plot(Range.rm,log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1))))
% hold on
% plot(Range.rm,log(Model.transmission(:,p_point(1)).^2))
% 
% %TT = Stemp(minJJ) + lapse(minII)/1000 .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
% %PP = Spress(minKK) .* (Stemp(minJJ)./TT).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
% %alpha = absorption_O2_770_model(TT,PP,Spectrum.nu_online(1),Model.WV(:,p_point)); %[m-1] Funcrtion to calculate theoretical absorption
% %Transmission = exp(-cumtrapz(Range.rm,alpha));
% plot(Range.rm,log(Transmission(:,p_point(1)).^2)+const(p_point(1)))
% hold off
% legend('Counts','model','fit')
% 
% figure
% plot(Range.rm,exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Transmission(:,p_point(1)).^2)-const(p_point(1))))
% hold on
% plot(Range.rm,exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Model.transmission(:,p_point(1)).^2)))
% %plot(Range.rm,smooth(exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Transmission(:,p_point(1)).^2)-const(p_point(1))),5))
% plot(Range.rm,overlapCorr(:,p_point))
% hold off
% legend('Fit','model')


figure(3567)
plot(diag(alpha_total_raw(:,p_point)-alpha_0(:,p_point)),rkm)
xlim([-2e-4 2e-4])
ylim([0 5])
grid on
title(sprintf(['alpha final - alpha 0\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')






%=bg
if ishandle(6488)
    figure(6488)
    clf(6488)
end
figure(6488)
line(thr,Counts.bg_o2on-Counts.bg_o2off,'Color','r')
line(thr,Counts.bg_o2on_mol-Counts.bg_o2off_mol,'Color','g')
line(thr,Counts.bg_o2on_mol-Counts.bg_o2on,'Color','g','LineStyle','--')
line(thr,Counts.bg_o2off_mol-Counts.bg_o2off,'Color','r','LineStyle','--')
ylim([-5 5]);
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('Online Background - offline background')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(thr,T_finalm(9,:),'Parent',ax2,'Color','k')
line(thr,T_final_test(13,:),'Parent',ax2,'Color','b')
ylabel('Temperature horizontal')

if ishandle(64888)
    figure(64888)
    clf(64888)
end
figure(64888)
line(thr,smooth(-Counts.bg_o2on+Counts.bg_o2off,10),'Color','r')
ylim([-5 5])
%ylim([-25 25])
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('-Online Background + offline background')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(thr,smooth(alpha_0(14,:),10),'Parent',ax2,'Color','k')
line(thr,smooth(alpha_0(9,:),10),'Parent',ax2,'Color','b')
ylim([0 5e-4])
grid on
ylabel('alpha 0')

figure(3927)
plot(Data.Power.O2Online.TimeStamp,Data.Power.O2Online.LaserPower)
hold on
plot(Data.Power.O2Offline.TimeStamp,Data.Power.O2Offline.LaserPower)
%plot(Data.Power.O2Online.TimeStamp,smooth(Data.Power.O2Online.LaserPower,1000))
%plot(Data.Power.O2Offline.TimeStamp,smooth(Data.Power.O2Offline.LaserPower,1000))
xline(Time.ts(p_point(1))/60/60)
hold off
xlim([Time.ts(1) Time.ts(end)]/60/60)
grid on
title('Power Measuremet')
legend('Online','Offline')

figure
Model.T_pot = Model.T.*(0.9869232667160128./Model.P).^(0.286);
T_pot = T_finalm.*(0.9869232667160128./Model.P).^(0.286);
plot(diag(Model.T_pot(:,p_point)),rkm)
hold on
plot(diag(T_pot(:,p_point)),rkm)
hold off
legend('model','measurement')
grid on
xlabel('Potential Temperature (K)')
ylabel('Range (rkm)')

figure
Model.T_pot = Model.T.*(0.9869232667160128./Model.P).^(0.286);
T_pot = T_finalm.*(0.9869232667160128./Model.P).^(0.286);
plot(diff(diag(Model.T_pot(:,p_point)))/Range.rangeBin ,rkm(1:end-1))
hold on
plot(diff(diag(T_pot(:,p_point)))/Range.rangeBin,rkm(1:end-1))
xlim([-.15 .15])
hold off
legend('model','measurement')
grid on
xlabel('Potential Temperature lapse rate (K/m)')
ylabel('Range (rkm)')

