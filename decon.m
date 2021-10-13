clear all
clear all
x = zeros(50,1);
x(10:15)=1.5;
%x(1:6)=1.5;

%y = zeros(25,1);
%y(1:2)=2;
y=[1;1;2; 2; 2; 1];

%xy = conv(x,y,'same');
xy = conv(x,y,'valid');
% ftxy = fft(xy);
% y = padarray(y,(length(xy)-length(y))/2);
% fty = fft(y);
% 
% 
% xyDecon = ifft(ftxy./fty);
%xyDecon = ifftshift(ftxy./fty);

figure
plot(x)
hold on
plot(y)
plot(xy)
%plot(xyDecon,'.-')
xyDecon = fdeconv(xy,y);
[xyDecon2,xyDecon2R] = deconv(xy,y);
plot(xyDecon,'*-')
plot(xyDecon,'o-')
hold off
legend('x','y','xy','xyDecon')
