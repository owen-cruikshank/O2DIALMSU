function [N_Corr,errorEst] = countDeconvolution(N_pulse,N_model,rm,pulseLength)
%%countDeconvolution
%Deconvolve counts from a rectangular pulse
%countDeconvolution.m
%author: Owen Cruikshank
%date: 8/27/2020

    d_r = rm(2)-rm(1);
    pulseBins = floor(pulseLength/d_r);
    
    dN_pulse = diff(N_pulse,1,1)/d_r;
    dN_pulse(end+1,:) = dN_pulse(end,:);
    
    %two point central difference
%    dN_pulse(1) = (N_pulse(2)-N_pulse(1))/d_r;
%    dN_pulse(2:length(N_pulse)-1)=(N_pulse(3:end)-N_pulse(1:end-2))/2/d_r;
%    dN_pulse(length(N_pulse))=(N_pulse(end)-N_pulse(end-1))/d_r;
    
%     rmCorr=0;
%     rmCorrInd=1;
%     % N_on_Corr(rmCorrInd,1:2) = 0;
%      N_off_Corr(rmCorrInd,1:2) = 0;
%    figure()
%    hold on
%    N_on_Corr(1:pulseBins)=N_on(1:pulseBins,1);
    %%%N_Corr(1:pulseBins,:)=ones(pulseBins,:);
    N_Corr(1:pulseBins,1:length(N_pulse(1,:)))=N_model(1:pulseBins,:).*ones(pulseBins,length(N_pulse(1,:)));
    %N_off_Corr(1:pulseBins)=0;
    %for j=pulseBins+1:length(rm)
     for j=pulseBins+1:length(rm)
       N_Corr(j,:)=0;
%       Q = floor(rm(j)/(pulseLength+3*d_r));
       Q = floor(rm(j)/(pulseLength));
       for i=0:Q
          % Q
          % j
          % i
          % pulseBins
           dNInd = j-i*pulseBins;
           N_Corr(j,:) = dN_pulse(j-i*pulseBins,:)+N_Corr(j,:);

%            if j==length(rm)-2*pulseBins-3
%            plot(i,dN_pulse(j-i*pulseBins,2),'o')
%            end
    %        i
    %        rmInd = rm(j-i*pulseBins)
       end
       %j
       %(Q+1)*pulseBins
       
       N_Corr(j,:)=N_Corr(j,:)*pulseLength+N_model(j-(Q)*pulseBins,:);
       %N_Corr(j)=N_Corr(j)*pulseLength;
      
    end
    
    %error estimateion
    %delta(z) = -(1/30)(dz)^4 P_actual IV
    
    errorEst = (-1/30)*(d_r)^4 * diff(N_pulse,4,1)/(d_r)^4;
end