% Calculating pertabative absorption of O2
% author: Owen Cruikshank
% date: 11/30/2020

function [alpha_1, alpha_2] = pertAbsorption(alpha_0, T_etalon, T, P, rm,ts,rkm,m_air, nu_online, nu_scan_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,t_avg,c,kb,altitude)
% --- Spectral distribution using the initial temperature profile guess ---

    cB = 1.2;%Brullouion correction to doppler gaussian half width

    cB = -0.01*(rkm+altitude) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km
    %cB = -0.01*(rkm+altitude) + 1.45;
    %cB = -0.01*(rkm+altitude) + 1.4;
    %cB=1;

    c_doppler_O2 = m_air*c^2./(8*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         
    %mo2 = 5.314E-26;            %[kg] Mass O2 molecule   
    %c_doppler_O2 = mo2*c^2./(8*pi*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
    % %  c_doppler_O2 = m_air*c^2./(8*pi*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
    %  doppler_O2_un_ret = ((c_doppler_O2./T).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    % Check if doppler_o2_ret is normalized to 1 when integrated across frequency
    doppler_o2_ret_check = trapz(doppler_O2_ret,3).*nuBin*100;              %[none]

    
    %%
    % Written by John Smith
% October 21st, 2010
% University of Colorado at Boulder, CIRES
% John.A.Smith@Colorado.EDU
% MATLAB version 7.10.0.59 (R2010a) 64-bit
% Adapted from "Coherent Rayleigh-Brillouin Scattering"
% by Xingguo Pan

% Sets and computes all relevant parameters
% for s6 and s7 models given in Xingguo Pan's
% dissertation entitled "Coherent Rayleigh-
% Brillouin Scattering", 2003, Princeton Univ.

% set the temperature of the gas
% % tem=292;
% % p_atm=3;
% % p_atm=0.6;
% % 
% % tem=290.9766;
% % p_atm=.9674;

% create the domain of xi
% % N=1000;
% % N=301;
% % xi_lef=-5.0d0;
% % xi_rgt=5.0d0;
% % xi=linspace(xi_lef,xi_rgt,N);
% *** 
% 'xi' is dimensionless and should be scaled by k*v0/(2*pi) [Hz],
% where 'v0' is the most probable gas velocity, or sqrt(2*kb*T/m),
% and 'k' is 4*pi/lambda*sin(scatter_angle/2)
% ***

% % % % disp('Calculating RB spectrum')
% % % % 
% % % % 
% % % % % set fundamental constants
% % % % kb=1.3806503e-23;
% % % % 
% % % % % set laser parameters
% % % % 
% % % % lambda =10^7/nu_online(1);
% % % % lambda = lambda*10^-9;
% % % % angle=(89.4)*(pi/180);
% % % % k=sin(angle)*4*pi/lambda;
% % % % 
% % % % % set N2 gas quantities
% % % % m_m=(1.66053886e-27)*28.013;
% % % % viscosity=17.63e-6;
% % % % bulk_vis=viscosity*0.73;
% % % % thermal_cond=25.2e-3;
% % % % c_int=1.0;
% % % % 
% % % % freq = nu_scan_3D_short * 100 * 3e8 - nu_online(1) *100 * 3e8;
% % % % 
% % % % c_tr=3/2;
% % % % gamma_int=c_int/(c_tr+c_int);
% % % % eukenf=m_m*thermal_cond/(viscosity*kb*(c_tr+c_int));
% % % % 
% % % % rlx_int=1.5*bulk_vis/(viscosity*gamma_int);
% % % % tic
% % % % for i = 1:i_range
% % % %    parfor j = 1:i_time
% % % %         
% % % % % for i = 1:1
% % % % %     for j = 1:1
% % % %         
% % % %         p_atm = P(i,j);
% % % %         tem = T(i,j)
% % % % 
% % % % 
% % % % 
% % % % % compute most probable gas velocity
% % % % v0=sqrt(2*kb*tem/m_m);
% % % % 
% % % % % convert pressures and densities
% % % % p_pa=p_atm*1.01325e5;
% % % % p_torr=p_pa*0.00750061683;
% % % % n0=p_pa/(tem*kb);
% % % % 
% % % % % compute and set RBS model input parameters
% % % % 
% % % % y=n0*kb*tem/(k*v0*viscosity);
% % % % 
% % % % xi = freq/((sin(angle)*4*pi/lambda)*v0/(2*pi));
% % % % 
% % % % xi = permute(xi,[2 3 1]);
% % % % 
% % % % % run the code
% % % % % %[cohsig6,sptsig6]=crbs6(y,rlx_int,eukenf,c_int,c_tr,xi);
% % % % % % [cohsig7,sptsig7]=crbs7(y,rlx_int,eukenf,c_int,c_tr,xi);
% % % % % % 
% % % % % % RBcoh(i,j,:)=permute(cohsig7,[3 2 1]);
% % % % % % RBspt(i,j,:)=permute(cohsig7,[3 2 1]);
% % % % 
% % % % 
% % % %     end
% % % % end
% % % % load('cohsig7.mat','RBcoh','RBspt')
% % % % 
% % % % 
% % % % norm_O2_ret = trapz(RBspt,3).*nuBin*100;                   %[none] Lineshape integral
% % % % %    doppler_O2_ret = RBspt./norm_O2_ret;                       %[m] Normalized doppler lineshape
% % % %     
% % % %     save('cohsig7.mat','RBcoh','RBspt','doppler_O2_ret','doppler_O2_un_ret')
% % % %     toc


%%%%%load('cohsig7.mat','RBcoh','RBspt','doppler_O2_ret','doppler_O2_un_ret')

% % freq = xi*(sin(angle)*4*pi/lambda)*v0/(2*pi);
% % lambda_s=freq*lambda^2/3e8 + lambda;
% % %.2558*
% % sptsig7_2=sptsig7/(trapz(sptsig7)*(xi(2)-xi(1)));
% % cohsig7_2=cohsig7/(trapz(cohsig7)*(xi(2)-xi(1)));
% % 
% % sptsig7_3=sptsig7/(trapz(sptsig7)*(lambda_s(2)-lambda_s(1)));
% % cohsig7_3=cohsig7/(trapz(cohsig7)*(lambda_s(2)-lambda_s(1)))

% % % % norm_O2_ret = trapz(RBcoh,3).*nuBin*100;                   %[none] Lineshape integral
% % % % doppler_O2_ret = RBcoh./norm_O2_ret;                       %[m] Normalized doppler lineshape
% % % % 
% % % % % norm_O2_ret = trapz(RBspt,3).*nuBin*100;                   %[none] Lineshape integral
% % % % % doppler_O2_ret = RBspt./norm_O2_ret;                       %[m] Normalized doppler lineshape
% % % % figure()
% % % % plot(permute(nu_scan_3D_short,[3 2 1]),permute(doppler_O2_ret(1,338,:),[3 2 1]))

    
    %%
    
    
    
    % --- Backscatter Lineshape g ---
    g1_m = 1./BSR .* doppler_O2_ret ;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
    g1_a = zeros(i_range,i_time,i_scan_3D_short);                       % Initalize aerosol lineshape
    for i = 1:i_time
        g1_a(:,i,online_index(i)) = (1 - 1./BSR(:,i))/ nuBin / 100 ; %[m] aerosol backscatter lineshape
    end
    g1 = g1_a + g1_m;                                                   %[m] Combined backscatter lineshape

    g1_check = trapz(g1,3).*nuBin*100;                                %[none] Check if integral of g1 is normalized to 1


    %derivative of lineshape dg/dr
    dg1_dr = (g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(rangeBin*oversample); %[none] Derivative over oversamped range
    %dg1_dr = (g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(rangeBin); %[none] Derivative over oversamped range
    dg1_dr = interp1(rm(ind_r_lo),dg1_dr,rm,'nearest',nan);         %[none] Make dg/dr the same size as g
    
    disp('SG derivative')
    tic
    [~,g] = sgolay(2,oversample+1);
    parfor j=1:i_time
        for k=1:i_scan_3D_short
            dg1_dr(:,j,k) = conv(g1(:,j,k), factorial(1)/(-rangeBin)^1 * g(:,2), 'same');
        end
    end
    toc
    
% % %     %derivative of lineshape dg/dr
% % %     for i=1:i_range-2*oversample
% % %         dg1_dr(i,:,:) = (-3*g1(i,:,:)+4*g1(i+oversample,:,:)-g1(i+2*oversample,:,:))./(2*rangeBin*oversample);
% % %     end
% % %     for i=i_range-2*oversample+1:i_range-oversample
% % %         dg1_dr(i,:,:) = (g1(i+oversample,:,:)-g1(i,:,:))./(rangeBin*oversample);
% % %     end
% % %     dg1_dr = interp1(rm(ind_r_lo),dg1_dr,rm,'nearest',nan);         %[none] Make dg/dr the same size as g
    
    %%%%dg1_dr = interp1(rm(ind_r_lo)+(rangeBin*oversample)/2,dg1_dr,rm,'nearest',nan);         %[none] Make dg/dr the same size as g

% %     dg1_dr = (g1(2:end,:,:) - g1(1:end-1,:,:)) ./(rangeBin); %[none] Derivative over oversamped range
% %     dg1_dr = interp1(rm(1:end-1),dg1_dr,rm,'nearest','extrap');         %[none] Make dg/dr the same size as g

figure()
plot(trapz(dg1_dr(:,338,:),3)*nuBin*100,rm)


    %%
    tic
    %[absorption_f,~,f] = absorption_O2_770_model_wavenumber(T,P,nu_scan_3D_short,WV); %[m] lineshape function 
    %[absorption_f] = absorption_O2_770_model_wavenumber(T,P,nu_scan_3D_short,WV); %[m] lineshape function 


    for i = 1:i_range
        parfor j=1:i_time
    absorption_f(i,j,:) = absorption_O2_770_model_wavenumber(T(i,j),P(i,j),nu_scan_3D_short(1,1,:),WV(i,j)); %[m-1] Funcrtion to calculate theoretical absorption
        end
    end
    %absorption_old = absorption_O2_770_model_wavenumber_old(T,P,nu_scan_3D_short);

    toc
    %%
    f = nan(size(absorption_f));
    for i = 1:i_time
        f(:,i,:) = absorption_f(:,i,:) ./ absorption_f(:,i,online_index(i));
        %f(:,i,:) = absorption_f(:,i,:) ./ max(absorption_f(:,i,:),[],3);    %[none] Normalize cross section to line center
        %f_old(:,i,:) = absorption_old(:,i,:) ./ absorption_old(:,i,online_index(i));
    end

     %%    
    % --- Zeroth Order Transmission ---
    Tm0 = exp(-cumtrapz(rm,alpha_0.*f,1));      %[none] Zeroth order transmission
    %Tm0 = exp(-cumtrapz(rm,alpha_0.*f,1));      %[none] Zeroth order transmission
    
    % ------ fitting alpha zero ---------
% % %     polyf = fittype('poly1');
% % %     exclusion = zeros(size(alpha_0(26:106,:)));
% % %     alpha_0_fit = nan(size(alpha_0));
% % % for i = 1:i_time
% % %     exclusion(:,i) = isnan(alpha_0(26:106,i));
% % %     if isnan(alpha_0(26,i))
% % %         
% % %     else   
% % %         alpha_0_fit_obj = fit(rm(26:106,1),alpha_0(26:106,i),polyf,'Exclude',exclusion(:,i));
% % %         alpha_0_fit(:,i) = alpha_0_fit_obj(rm);
% % %     end
% % % end

    % Integrand terms
    % Online
    zeta = g1.*T_etalon;                        %[m]
    eta = dg1_dr.*T_etalon;                     %[none]

    % Integrated terms
    % Online
    zeta_int = trapz(zeta.*Tm0,3)*nuBin*100;              %[none]
    eta_int = trapz(eta.*Tm0,3)*nuBin*100;                %[1/m]
    zeta_ls_int = trapz(zeta.*Tm0.*(1-f),3)*nuBin*100;    %[none]
    % Offline
    zeta2_int = trapz(zeta,3)*nuBin*100;                  %[none]
    eta2_int = trapz(eta,3)*nuBin*100;                    %[1/m]

    % === First Order ===
    W1 = zeta_ls_int./zeta_int;                 %[none]
    G1 = eta_int./zeta_int - eta2_int./zeta2_int;%[1/m]

    %alpha_1_raw = 0.5.*(alpha_0_fit.*W1 + G1);  %[1/m]
    alpha_1_raw = 0.5.*(alpha_0.*W1 + G1);      %[1/m]

    % Moving average
    
   k = ones(oversample,t_avg)./(oversample*t_avg);
   t_step = ts(2)-ts(1);
   rangeBin = rm(2)-rm(1);
%     alpha_1_pad = padarray(alpha_1_raw,[oversample/2,t_avg/2],'replicate');
%     alpha_1_filt = filter2(k,alpha_1_pad,'valid');
%     alpha_1 = interp2(ts-t_step/2,rm-rangeBin/2,alpha_1_filt(1:end-1,1:end-1),ts,rm);
%     alpha_1 = fillmissing(alpha_1,'nearest',1); % Fill in NaNs in dimension 1
%     alpha_1 = fillmissing(alpha_1,'nearest',2); % Fill in NaNs in dimension 2
   alpha_1 = alpha_1_raw;

    alpha_1s = zeros(i_range,i_time);
    % ----- smoothing alpha 1 -----
    for i = 1:i_time
        alpha_1s(:,i) = 1.15.*smooth(alpha_1(:,i), 4*oversample+1);
    end

    % --- First Order Transmission Tm1 ---
    Tm1 = exp(-cumtrapz(rm,oversample.*alpha_1.*f,1));      %[none] First order transmission
    %Tm1 = exp(-cumtrapz(rm,4*alpha_1.*f,1));      %[none] First order transmission

    % === Second Order ===
    % Integrated terms
    zeta_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1),3)*nuBin*100;             %[none]
    eta_Tm1_int = trapz(eta.*Tm0.*(1-Tm1),3)*nuBin*100;               %[1/m]
    zeta_ls_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1).*(1-f),3)*nuBin*100;   %[none]

    W2 = (zeta_ls_int.*zeta_ls_Tm1_int./(zeta_int.^2)) - (zeta_ls_Tm1_int./zeta_int);   %[none]
    G2 = (eta_int.*zeta_Tm1_int./(zeta_int.^2)) - (eta_Tm1_int./zeta_int);              %[1/m]

    %alpha_2_raw = 0.5.*(alpha_1.*W1 + alpha_0_fit.*W2 + G2);
    alpha_2_raw = 0.5.*(alpha_1.*W1 + alpha_0.*W2 + G2);    %[1/m]

    %Moving average
%     alpha_2_pad = padarray(alpha_2_raw,[oversample/2,t_avg/2],'replicate');
%     alpha_2_filt = filter2(k,alpha_2_pad,'valid');
%     alpha_2 = interp2(ts-t_step/2,rm-rangeBin/2,alpha_2_filt(1:end-1,1:end-1),ts,rm);
%     alpha_2 = fillmissing(alpha_2,'nearest',1); % Fill in NaNs in dimension 1
%     alpha_2 = fillmissing(alpha_2,'nearest',2); % Fill in NaNs in dimension 2
   alpha_2=alpha_2_raw;
%%
    figure()
    plot(G1(:,338),rm)
    title('G1')
    
    figure()
    plot(G2(:,338),rm)
    title('G2')
    
        figure()
    plot(W1(:,338),rm)
    title('W1')
    
    figure()
    plot(W2(:,338),rm)
    title('W2')
%%

end