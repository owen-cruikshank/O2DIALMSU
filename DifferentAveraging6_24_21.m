
path = 'C:\Users\oencr\OneDrive - Montana State University\research s19\Reports\6_24_21\';
path = 'C:\Users\Owen C\OneDrive - Montana State University\research s19\Reports\6_24_21\';
load([path 'int10min70m_Aavg60min300m_Tavg60min300m.mat'])
T_finalm1 = T_finalm;
load([path 'int10min70m_Aavg40min300m_Tavg40min300m.mat'])
T_finalm2 = T_finalm;
load([path 'int10min37m_Aavg40min300m_Tavg40min300m.mat'])
T_finalm3 = T_finalm;
T_sonde2=Sonde.T_sonde;
rm2 = Range.rm;
load([path 'int10min112m_Aavg60min300m_Tavg60min300m.mat'])
T_finalm5 = T_finalm;
rm3 = Range.rm;
load([path 'int10min150m_Aavg40min300m_Tavg40min300m.mat'])
T_finalm6 = T_finalm;
rm6 = Range.rm;
load([path 'int10min70m_Aavg40min300m_Tavg40min300mCorrected.mat'])
T_finalm4 = T_finalm;

%%

sonde_index = 16;
p_point = Sonde.sonde_ind(:,sonde_index);

figure
plot(T_sonde2(:,sonde_index),rm2)
hold on
plot(T_finalm1(:,p_point(1)),Range.rm,'--')
plot(T_finalm2(:,p_point(1)),Range.rm,'.-')
plot(T_finalm3(:,p_point(1)),rm2,'--')
plot(T_finalm4(:,p_point(1)),Range.rm,'linewidth',2)
plot(T_finalm5(:,p_point(1)),rm3)
plot(T_finalm6(:,p_point(1)),rm6,'--','linewidth',2)
hold off
legend('Sonde','int10min70m_Aavg60min300m_Tavg60min300m','int10min70m_Aavg40min300m_Tavg40min300m','int10min37m_Aavg40min300m_Tavg40min300m','int10min70m_Aavg40min300m_Tavg40min300mCorrected')