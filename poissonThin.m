function Counts = poissonThin(Counts)

disp('Thining profiles')
tic
%====Thin counts into two profiles====
%fon = binornd(round(Counts.o2on.*Counts.NBins),0.5);
fon = TEST_SUB_mybinornd( round((Counts.o2on+Counts.bg_o2on).*Counts.NBins), 0.5);
gon = round((Counts.o2on+Counts.bg_o2on).*Counts.NBins)-fon;
%foff = binornd(round(Counts.o2off.*Counts.NBins),0.5);
foff = TEST_SUB_mybinornd( round((Counts.o2off+Counts.bg_o2off).*Counts.NBins), 0.5);
goff = round((Counts.o2off+Counts.bg_o2off).*Counts.NBins)-foff;
%fon_mol = binornd(round(Counts.o2on_mol.*Counts.NBins),0.5);
fon_mol = TEST_SUB_mybinornd( round((Counts.o2on_mol+Counts.bg_o2on_mol).*Counts.NBins), 0.5);
gon_mol = round((Counts.o2on_mol+Counts.bg_o2on_mol).*Counts.NBins)-fon_mol;
%foff_mol = binornd(round(Counts.o2off_mol.*Counts.NBins),0.5);
foff_mol = TEST_SUB_mybinornd( round((Counts.o2off_mol+Counts.bg_o2off_mol).*Counts.NBins), 0.5);
goff_mol = round((Counts.o2off_mol+Counts.bg_o2off_mol).*Counts.NBins)-foff_mol;
toc

%=====Find background of thinned profiles=====
fon_bg = (Counts.bg_o2on.*Counts.NBins)/2;% Take mean of last data points
fon_mol_bg = (Counts.bg_o2on_mol.*Counts.NBins)/2;% Take mean of last data points
foff_bg = (Counts.bg_o2off.*Counts.NBins)/2;% Take mean of last data points
foff_mol_bg = (Counts.bg_o2off_mol.*Counts.NBins)/2;% Take mean of last data points

%=====Find optimal filters====
[Ezon,Eton,rangeWidthon,timeWidthon] = findMinE(fon,gon,fon_bg);
[Ezoff,Etoff,rangeWidthoff,timeWidthoff] = findMinE(foff,goff,foff_bg);
[Ezon_mol,Eton_mol,rangeWidthon_mol,timeWidthon_mol] = findMinE(fon_mol,gon_mol,fon_mol_bg);
[Ezoff_mol,Etoff_mol,rangeWidthoff_mol,timeWidthoff_mol] = findMinE(foff_mol,goff_mol,foff_mol_bg);

Counts.Poissonthin.timeWidthon = timeWidthon;
Counts.Poissonthin.timeWidthoff = timeWidthoff;
Counts.Poissonthin.timeWidthon_mol = timeWidthon_mol;
Counts.Poissonthin.timeWidthoff_mol = timeWidthoff_mol;
Counts.Poissonthin.rangeWidthon = rangeWidthon;
Counts.Poissonthin.rangeWidthoff = rangeWidthoff;
Counts.Poissonthin.rangeWidthon_mol = rangeWidthon_mol;
Counts.Poissonthin.rangeWidthoff_mol = rangeWidthoff_mol;

Counts.Poissonthin.Ezon=Ezon;
Counts.Poissonthin.Eton=Eton;
Counts.Poissonthin.Ezoff=Ezoff;
Counts.Poissonthin.Etoff=Etoff;
Counts.Poissonthin.Ezon_mol=Ezon_mol;
Counts.Poissonthin.Eton_mol=Eton_mol;
Counts.Poissonthin.Ezoff_mol=Ezoff_mol;
Counts.Poissonthin.Etoff_mol=Etoff_mol;

%==== Appy optimal filters to Count data====
%==Filter in range==
for iii = 1:size(Counts.o2on,2)    
        nz = round(4*rangeWidthon(iii)); %number of grid points 
        z = (-nz:nz)';%filter grid
        %kern = gaussmf(z,[rangeWidthon(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/rangeWidthon(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2on(:,iii))),kern,'same');
        %Counts.o2on(:,iii) = conv2(Counts.o2on(:,iii)-Counts.bg_o2off(iii),kern,'same')./norm;
        %Counts.o2on(:,iii) = conv2(Counts.o2on(:,iii),kern,'same')./norm;
        Counts.o2on(:,iii) = nanconv(Counts.o2on(:,iii),kern,'edge','nanout');
        
        nz = round(4*rangeWidthoff(iii)); %number of grid points 
        z = (-nz:nz)';%filter grid
        %kern = gaussmf(z,[rangeWidthoff(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/rangeWidthoff(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2off(:,iii))),kern,'same');
        %Counts.o2off(:,iii) = conv2(Counts.o2off(:,iii)-Counts.bg_o2off(iii),kern,'same')./norm;
        %Counts.o2off(:,iii) = conv2(Counts.o2off(:,iii),kern,'same')./norm;
        Counts.o2off(:,iii) = nanconv(Counts.o2off(:,iii),kern,'edge','nanout');
        
        nz = round(4*rangeWidthon_mol(iii)); %number of grid points 
        z = (-nz:nz)';%filter grid
        %kern = gaussmf(z,[rangeWidthon_mol(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/rangeWidthon_mol(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2on_mol(:,iii))),kern,'same');
        %Counts.o2on_mol(:,iii) = conv2(Counts.o2on_mol(:,iii)-Counts.bg_o2off_mol(iii),kern,'same')./norm;
        %Counts.o2on_mol(:,iii) = conv2(Counts.o2on_mol(:,iii),kern,'same')./norm;
        Counts.o2on_mol(:,iii) = nanconv(Counts.o2on_mol(:,iii),kern,'edge','nanout');
        
        nz = round(4*rangeWidthoff_mol(iii)); %number of grid points 
        z = (-nz:nz)';%filter grid
        %kern = gaussmf(z,[rangeWidthoff_mol(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/rangeWidthoff_mol(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2off_mol(:,iii))),kern,'same');
        %Counts.o2off_mol(:,iii) = conv2(Counts.o2off_mol(:,iii)-Counts.bg_o2off_mol(iii),kern,'same')./norm; 
        %Counts.o2off_mol(:,iii) = conv2(Counts.o2off_mol(:,iii),kern,'same')./norm; 
        Counts.o2off_mol(:,iii) = nanconv(Counts.o2off_mol(:,iii),kern,'edge','nanout');
end
%==Filter in time==
for iii = 1:size(Counts.o2on,1)    
        nz = round(4*timeWidthon(iii)); %number of grid points 
        z = (-nz:nz);%filter grid
        %kern = gaussmf(z,[timeWidthon(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/timeWidthon(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2on(iii,:))),kern,'same');
        %Counts.o2on(iii,:) = conv2(Counts.o2on(iii,:),kern,'same')./norm;
        Counts.o2on(:,iii) = nanconv(Counts.o2on(:,iii),kern,'edge','nanout');
        
        nz = round(4*timeWidthoff(iii)); %number of grid points 
        z = (-nz:nz);%filter grid
        %kern = gaussmf(z,[timeWidthoff(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/timeWidthoff(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2off(iii,:))),kern,'same');
        %Counts.o2off(iii,:) = conv2(Counts.o2off(iii,:),kern,'same')./norm;
        Counts.o2off(:,iii) = nanconv(Counts.o2off(:,iii),kern,'edge','nanout');
        
        nz = round(4*timeWidthon_mol(iii)); %number of grid points 
        z = (-nz:nz);%filter grid
        %kern = gaussmf(z,[timeWidthon_mol(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/timeWidthon_mol(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2on_mol(iii,:))),kern,'same');
        %Counts.o2on_mol(iii,:) = conv2(Counts.o2on_mol(iii,:),kern,'same')./norm;
        Counts.o2on_mol(:,iii) = nanconv(Counts.o2on_mol(:,iii),kern,'edge','nanout');
        
        nz = round(4*timeWidthoff_mol(iii)); %number of grid points 
        z = (-nz:nz);%filter grid
        %kern = gaussmf(z,[timeWidthoff_mol(iii),0]);%fitler kernel in range
        kern = exp(-z.^2/timeWidthoff_mol(iii).^2);
        kern = kern/sum(sum(kern));
        %norm = conv2(ones(size(Counts.o2off_mol(iii,:))),kern,'same');
        %Counts.o2off_mol(iii,:) = conv2(Counts.o2off_mol(iii,:),kern,'same')./norm; 
        Counts.o2off_mol(:,iii) = nanconv(Counts.o2off_mol(:,iii),kern,'edge','nanout');
end
end

function [Ez,Et,minSigz,minSigt] = findMinE(f,g,bg)
    disp('creating filter')
    %====find best filter in range====
    filt_size = logspace(-1,1.5,40); %create filter size in terms of grid points
    Ez = ones(1,length(f(1,:)),length(filt_size));
    for jj = 1:length(filt_size)
        fprintf('Range %f',jj)
        nz = round(4*filt_size(jj)); %number of grid points 
        z = (-nz:nz)';%filter grid
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
        Ez(:,:,jj) = sum(fFilt+bg-g.*log(fFilt+bg),1,'omitnan');
    end
    [~,minEind]=min(Ez,[],3);
    minSigz = ones(1,size(f,2));
    for ii = 1:size(f,2)
        minSigz(ii) = filt_size(minEind(ii));
    end

    %====find best filter in time====
    filt_size = logspace(-1.5,5,100);
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

function [ res ] = TEST_SUB_mybinornd( N, p )
    [row_cnt, col_cnt] = size(N);
    res = zeros(row_cnt, col_cnt);
    for ii=1:row_cnt
       for jj=1:col_cnt
           if isnan(N(ii,jj))
               res(ii,jj)=NaN;
           else
               res(ii, jj) = sum(rand(1,N(ii,jj))<p);
           end
       end
    end
end
    
    