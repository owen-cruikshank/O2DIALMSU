% ray tracing analysis
clear all
% rt =.356/2;% telescope radius
% r0 = .0508;% beam radius at telescope
% theta = 5.3e-6;% beam divergence

rt =.406/2;% telescope radius
r0 = .180/2;% beam radius at telescope
theta = 70e-6/2;% beam divergence
%theta = 5.3e-6;% beam divergence
theta = 2*770e-9/pi/114e-3; %half angle divergence, M2*lambda/pi/effective beam diameter

R =0:50:12000;%range vector
rb = r0+R.*theta;

dt = [0 50 100 150 200 300 500 1000]*10^-6;

R1 = [zeros(size(R))
    rt./R ];
R2 = [rb
    (rt-rb)./R ];
R3 = [rb
    zeros(size(R))];
R4 = [rb
    -(rt+rb)./R ];

% s{1} = [1 R(ii)
%         0 1];
s{2} = [1 0
    -1/1.218 1];
s{3} = [1 1.218+.060
    0 1];
s{4} = [1 0
    -1/.060 1];
s{5} = [1 .1
       0 1];
s{6} = [1 0
    -1/.080 1];
s{7} = [1 .080-.02
    0 1];
s{8} = [1 0
    1/.02 1];
s{9} = [1 .02
    0 1];
s{10} = [1 0
    -1/.011 1];
% s{11} = [1 .011+dt(jj)
%         0 1];
s{11} = [1 .011+dt(1)
        0 1];


for jj = 1:length(R)
    RR1{1}(:,jj) = R1(:,jj); 
    RR2{1}(:,jj) = R2(:,jj);
    RR3{1}(:,jj) = R3(:,jj);
    RR4{1}(:,jj) = R4(:,jj);
    s{1} = [1 R(jj)
            0 1];
    ind = 1;
    for ii = 1:length(s)
        RR1{ii+1}(:,jj) = s{ii}*RR1{ii}(:,jj);
        RR2{ii+1}(:,jj) = s{ii}*RR2{ii}(:,jj);
        RR3{ii+1}(:,jj) = s{ii}*RR3{ii}(:,jj);
        RR4{ii+1}(:,jj) = s{ii}*RR4{ii}(:,jj);

        
        if rem(ii, 2) == 1 %is odd
            d(ind,jj) = s{ii}(1,2);
            h1(ind,jj) = RR1{ii+1}(1,jj);
            h2(ind,jj) = RR2{ii+1}(1,jj);
            h3(ind,jj) = RR3{ii+1}(1,jj);
            h4(ind,jj) = RR4{ii+1}(1,jj);
            ind=ind+1;
        end
    end
%     d(ind+1,jj) = s{ii}(1,2);
%     h1(ind+1,jj) = RR1{ii+1}(1,jj);
%     h2(ind+1,jj) = RR2{ii+1}(1,jj);
%     h3(ind+1,jj) = RR3{ii+1}(1,jj);
%     h4(ind+1,jj) = RR4{ii+1}(1,jj);
end


lensRadius(1) = 406/2;
lensRadius(2) = 25.4/2;
lensRadius(3) = 25.4/2;
lensRadius(4) = 12.7/2;
lensRadius(5) = 4.4/2;
lensRadius(6) = .105/2;
lensRadius = lensRadius/1000;

x = cumsum(d,1) -d(1,:);
%%
figure
plot(x(:,end),h1(:,end),'O-')
hold on
plot(x(:,end),h2(:,end),'O-')
plot(x(:,end),h3(:,end),'O-')
plot(x(:,end),h4(:,end),'O-')
plot(x(:,end),lensRadius,'*')
plot(x(:,end),-lensRadius,'*')

plot(x(end-1,end),5.5/2/1000,'*')
plot(x(end-1,end),-5.5/2/1000,'*')
legend('1','2','3','4')
hold off
grid on

