syms P T nu0 dair S0 gL ep nair nu t wv

t = -10:.2:10;

T0 = 296;
P0 = 1;

c = 2.99792458E8;                       %[m/s] speed of light 
kB = 1.38065E-23;                       %[J/K] Boltzman's constant 
h = 6.62607004E-34;                     %[Js] Planck's constant
q_O2 = .2095;                           %[unitless] O2 atmospheric mixing ratio 
mo2 = 5.313525282632483E-26;            % Mass O2 molecule [kg]

nu_s = nu0 +dair*P;

ST = S0*(T0/T)*exp((h*c)/kB*(1/T0-1/T)*ep);
gLT = gL*(P/P0)*(T0/T)^nair;
gDT = nu_s*sqrt(2*kB*T*log(2)/mo2)/c;
x = (nu-nu_s)/gDT *sqrt(log(2));
y = (gLT/gDT)*sqrt(log(2));
k = ST/gDT *sqrt(log(2)/pi);

%%V = y/pi * int(exp(-t.^2)./(y^2+(x-t).^2),t);

V = y/pi * trapz(t,exp(-t.^2)./(y^2+(x-t).^2));

alpha = V*k*(P/(kB*T)-wv)*q_O2;

alphaDT = diff(alpha,T);
alphaDP = diff(alpha,P);
alphaDWV = diff(alpha,wv);

dalpha = sqrt((dp*alphaDP).^2+(dt*alphaDP).^2)
