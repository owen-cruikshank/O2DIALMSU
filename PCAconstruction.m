%principal component analysis construction
clear all
N=100;
%N=1000;
%T vector
Tmin = 140;
Tmax = 320;
T = linspace(Tmin,Tmax,N);%[K]
%P vector
Pmin = .1;
Pmax = 1.2;
P = linspace(Pmin,Pmax,N);%[atm]

%freqency
lambdaCenter = 769.7958;
lambdaWidth = 0.2;
lambdaDelta = .00002;
lambda = lambdaCenter-lambdaWidth/2:lambdaDelta:lambdaCenter+lambdaWidth/2;%[nm]
nu = 10^7./lambda;


%%
lambda_online = 769.7958;
lambda_offline = 770.1085;

%lambda_online(1:length(ts)) = 769.7954;
%lambda_offline(1:length(ts)) = 770.1081;
nu_online = 10^7./lambda_online;                    %[cm-1] Online wavenumber
nu_offline = 10^7./lambda_offline;                  %[cm-1] Offline wavenumber

nu01 = nu_online;                                   %[cm-1] Set center of scan to 
nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
nuBin = 0.00222;                                    %[cm-1] Scan increment
%nuBin = nuBin/2;
%nuBin = 8.8800e-04;
%nuBin = 0.0007;                                    %[cm-1] Scan increment
nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nu01_off = nu_offline;
nuMin_off = nu01_off-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = nu01_off+0.334;                                 %[cm-1] Scan upper bound
nu_scan_off = (nuMin_off:nuBin:nuMax_off);

% %new nu
% lambdaCenter = 769.7958;
% lambdaWidth = 0.2;
% lambdaDelta = .00002;
% lambda = lambdaCenter-lambdaWidth/2:lambdaDelta:lambdaCenter+lambdaWidth/2;%[nm]
% nu_scan = 10^7./lambda;
% %nu_scan = permute(nu,[3 2 1]);


%i_scan = length(nu_scan);                           %[none] length of scan vector

lambda_scan = 10^7./nu_scan;                        %[nm](1 x lambda) Scan vector
%f_scan = nu_scan * c * 100;                         %[Hz]( x f) Scan vector

nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension

lambda_scan_3D_short = 10^7./nu_scan_3D_short;
lambda_scan_3D_short_off = 10^7./nu_scan_3D_short_off;
i_scan_3D_short = length(nu_scan_3D_short);         %[none] length of scan vector

del_nu = nu_scan_3D_short-nu_online;                %[1/cm] difference from center
del_lambda = lambda_scan_3D_short-lambda_online;
nu = nu_scan_3D_short;
lambda = lambda_scan;

% nu = permute(nu01,[3 2 1]);
% lambda = permute(lambda_online,[3 2 1]);
%%


% nu01 = 1.299045834338435e+04;                                   %[cm-1] Set center of scan to 
% nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
% nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
% nuBin = 0.00222;                                    %[cm-1] Scan increment
% %nuBin = 0.001;                                    %[cm-1] Scan increment
% nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
% %lambda = 10^7./nu_scan;
% %i_scan = length(nu_scan);                           %[none] length of scan vector
% 
% %nu = permute(nu_scan,[3 1 2]);
% nu = permute(nu,[3 1 2]);
% nuCenter = 10^7./lambdaCenter;


%construct P and T as j
TP = zeros(N*N,2);
increment = 1;
for i = 1:N
    for j = 1:N
        TP(increment,1) = T(i);
        TP(increment,2) = P(j);
        increment = increment+1;
    end
end
muT = mean(T);
sigmaT = std(T);
muP = mean(P);
sigmaP = std(P);

%%
%create training matrix
%y = zeros(1,N*N,length(lambda));
vw = zeros(N*N,1,1);
q = 10;
parfor i = 1:length(nu)
    %for j = 1:N*N
        [~,y(i,:),~] = absorption_O2_770_model(TP(:,1),TP(:,2),nu(i),0);
        %[a,y,~] = absorption_O2_770_model_wavenumber(TP(:,1),TP(:,2),nu,0);

    %end
    %y = permute(y,[3 2 1]);
end
%%
clear theta
No = 20;

for j = 1:N*N
    normT = (TP(j,1)-muT)/sigmaT;
    normP = (TP(j,2)-muP)/sigmaP;
    inc=1;
    for n = 1:No
        for m = 1:n
            theta(j,inc) = normT^(n-m) * normP^(m-1);
            inc=inc+1;
        end
    end
end
theta = theta';

%%
clear muY
muY = mean(y,2);

%%
clear U S V
[U,S,V] = svd(y);

%%
clear W
W = U'*(y-muY);

%%
clear C
C = W*pinv(theta);
%%
clear M
M = U*C;

%%
%calculating new
clear ynew
for j = 1:N*N
    ynew(:,j) = muY + M*theta(:,j);
end

%%
%save
%save('PCA_1_11_21singleOnline.mat','muY','M','muT','muP','sigmaT','sigmaP','nu')
save('PCA_1_11_21single.mat','muY','M','muT','muP','sigmaT','sigmaP','nu')
%%

%load('PCA_5_18_20.mat')

%%
error = abs((y-ynew)./y);

%%
TPind = 7583;
figure()
plot(lambda,y(:,TPind));
hold on
plot(lambda,ynew(:,TPind))
legend('y','ynew')
hold off

figure()
plot(lambda,error(:,TPind-5:TPind+5))

% figure(3)
% imagesc(TP(:,1),TP(:,2),error)

%%
%arbitrary
Ts = 1.2 + 273.15;          %surface temperature from weather station [K]
Ps = 839 / 1013.25;         %absolute surface pressure from weather station [atm]
tic
parfor i = 1:length(nu)
    %for j = 1:N*N
        [~,yabs(i,:),~] = absorption_O2_770_model(Ts,Ps,nu(i),0);
        %[a,y,~] = absorption_O2_770_model_wavenumber(TP(:,1),TP(:,2),nu,0);

    %end
    %y = permute(y,[3 2 1]);
end
toc
clear thetaNew
tic

    normT = (Ts-muT)/sigmaT;
    normP = (Ps-muP)/sigmaP;
    inc=1;
    No = 20;
    for n = 1:No
        for m = 1:n
            thetaNew(1,inc) = normT^(n-m) * normP^(m-1);
            inc=inc+1;
        end
    end
    thetaNew = thetaNew';
ynewtest = muY + M*thetaNew;
toc


errornew = abs((yabs(1:end-9)-ynewtest(10:end))./yabs(1:end-9));

errornew = abs((yabs-ynewtest)./yabs);

figure()
plot(lambda,yabs);
hold on
plot(lambda,ynewtest)
legend('y','ynew')
hold off

figure()
plot(lambda,errornew*100)
ylabel('error %')



