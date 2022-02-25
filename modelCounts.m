
function [Model] = modelCounts(Counts,Model,HSRL,Time,Range,Spectrum)
%model photon returns and overlap function
%7/15/21
%Owen Cruikshank

%Constants
h = 6.62607004e-34;
c = 2.9989e8;   % Speed of light m/s

%instrument parameters
A = 1;%area of telescope
eta_O = 1;%overlap
eta_D = .6*(1/18);%detector
eta_R = 7.4516e-22 * 1.4122*2.5*5;%reciever
E_pulse_on = 5 * 10^-6;%pulse energy (J)
E_pulse_off = 5 * 10^-6;%pulse energy (J)

%--Outgoing photons
pulse_rate = 7000;%(Hz)
avg_time = 2;%(sec)
transmitPhotons_on = E_pulse_on./(h*c./Spectrum.lambda_online*10^-9) * (pulse_rate/2) * avg_time*60 ;
%transmitPhotons_off = E_pulse_off./(h*c./lambda_offline*10^-9) * (pulse_rate/2) * avg_time*60 ;

%--transmission lidar ratio using HSRL ans model
Sa = 60;
Ta = exp(-cumtrapz(Range.rm,HSRL.Ba*Sa,1));
Sm = (8/3)*pi;
Tm = exp(-cumtrapz(Range.rm,HSRL.Bm*Sm,1));

%--Lidar equation
Fon  = Tm.^2.*Ta.^2.*Model.transmission.^2.*(HSRL.Bm+HSRL.Ba)./Range.rm.^2;
Foff = Tm.^2.*Ta.^2.*(HSRL.Bm+HSRL.Ba)./Range.rm.^2;
Foff_mol = Tm.^2.*Ta.^2.*(HSRL.Bm)./Range.rm.^2;


Model.N_on  = transmitPhotons_on  .* eta_R .* eta_D .* A .* Range.rangeBin .* Fon;
Model.N_off = transmitPhotons_on .* eta_R .* eta_D .* A .* Range.rangeBin .* Foff;

%--Pulse convolution
pulseLength=(1e-6)*c/2;%(m)
pulseBins = round(pulseLength/Range.rangeBin);
pulseCorr=1/117.6 * 1/1.3 * 1.2313 * 1/1.2205 *2.0213; %scaling factor

%pulse width modeling
Model.N_on_pulse = zeros(size(Model.N_on));
Model.N_off_pulse = zeros(size(Model.N_off));
for i = 1:length(Range.rm)
    if i<=pulseBins
        %1:i
%         N_on_pulse(i) = pulseCorr* dr*trapz(transmitPhotons_on  .* eta_R .* eta_D .* A .* dr .* Fon(1:i));
%         N_off_pulse(i) = pulseCorr* dr*trapz(transmitPhotons_on  .* eta_R .* eta_D .* A .* dr .*Foff(1:i));
        Model.N_on_pulse(i,:) = pulseCorr* trapz(Range.rangeBin,Model.N_on(1:i,:),1);
        Model.N_off_pulse(i,:) = pulseCorr* trapz(Range.rangeBin,Model.N_off(1:i,:),1);
    else
        %i-pulseBins+1:i
%         N_on_pulse(i) = pulseCorr* dr*trapz(transmitPhotons_on  .* eta_R .* eta_D .* A .* dr .*Fon(i-pulseBins+1:i));
%         N_off_pulse(i) = pulseCorr* dr*trapz(transmitPhotons_on  .* eta_R .* eta_D .* A .* dr .*Foff(i-pulseBins+1:i));
        Model.N_on_pulse(i,:) = pulseCorr* trapz(Range.rangeBin,Model.N_on(i-pulseBins+1:i,:),1);
        Model.N_off_pulse(i,:) = pulseCorr* trapz(Range.rangeBin,Model.N_off(i-pulseBins+1:i,:),1);
    end
end

%%
%--Calculate overlap
Model.OverlapOn = Counts.o2on./Model.N_on;
Model.OverlapOff = Counts.o2off./Model.N_off;
Model.OverlapOn_pulse = Counts.o2on./Model.N_on_pulse;
Model.OverlapOff_pulse = Counts.o2off./Model.N_off_pulse;