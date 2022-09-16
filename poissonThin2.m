function Counts = poissonThin2(Counts,cloud_SDm_above)
%File: poissonThin.m
%Date: 02/4/2022
%Author: Owen Cruikshank
%Inputs:
%   -Counts, Structure
%       -Counts.o2on (range x time): integrated and background subtracted
%       o2 online combined photon counts
%       -Counts.o2off (range x time): integrated and background subtracted
%       o2 offline combined photon counts
%       -Counts.o2on_mol (range x time): integrated and background subtracted
%       o2 online molecular photon counts
%       -Counts.o2off_mol (range x time): integrated and background subtracted
%       o2 offline molecular photon counts
%       -Counts.NBins (): number of raw data bins integrated into each bin
%       -Counts.bg_o2on (1 x time): background counts for o2on
%       -Counts.bg_o2off (1 x time): background counts for o2off
%       -Counts.bg_o2on_mol (1 x time): background counts for o2on_mol
%       -Counts.bg_o2off_mol (1 x time): background counts for o2off_mol
%
%Outputs:
%   -Counts, Structure
%       -Counts.o2on (range x time): optimally smoothed, integrated, and background subtracted
%       o2 online combined photon counts
%       -Counts.o2off (range x time): optimally smoothed, integrated, and background subtracted
%       o2 offline combined photon counts
%       -Counts.o2on_mol (range x time): optimally smoothed, integrated, and background subtracted
%       o2 online molecular photon counts
%       -Counts.o2off_mol (range x time): optimally smoothed, integrated, and background subtracted
%       o2 offline molecular photon counts
    %   -Counts.Poissonthin
    %         -Counts.Poissonthin.timeWidthon (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.timeWidthoff (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.timeWidthon_mol (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.timeWidthoff_mol (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.rangeWidthon (1 x time): Gaussian filter width as number
    %         of bins in range
    %         -Counts.Poissonthin.rangeWidthoff (1 x time): Gaussian filter width as number
    %         of bins in range
    %         -Counts.Poissonthin.rangeWidthon_mol (1 x time): Gaussian filter width as number
    %         of bins in range
    %         -Counts.Poissonthin.rangeWidthoff_mol (1 x time): Gaussian filter width as number
    %         of bins in range
    %         
    %         -Counts.Poissonthin.Ezon=Ezon;
    %         -Counts.Poissonthin.Eton=Eton;
    %         -Counts.Poissonthin.Ezoff=Ezoff;
    %         -Counts.Poissonthin.Etoff=Etoff;
    %         -Counts.Poissonthin.Ezon_mol=Ezon_mol;
    %         -Counts.Poissonthin.Eton_mol=Eton_mol;
    %         -Counts.Poissonthin.Ezoff_mol=Ezoff_mol;
    %         -Counts.Poissonthin.Etoff_mol=Etoff_mol;

disp('Thining profiles')

% Counts.o2on(cloud_SDm_above) = nan;
% Counts.o2off(cloud_SDm_above) = nan;
% Counts.o2on_mol(cloud_SDm_above) = nan;
% Counts.o2off_mol(cloud_SDm_above) = nan;
% Counts.wvon(cloud_SDm_above) = nan;
% Counts.wvoff(cloud_SDm_above) = nan;


% rng(0,'simdTwister') %faster random number generator
tic
% thin profiles into set f.
% input is rounded, background added, and multiplied by interation to be
% close to raw counts.
Counts.fon = TEST_SUB_mybinornd( round((Counts.o2on+Counts.bg_o2on).*Counts.NBins), 0.5)./Counts.NBins;

Counts.foff  = TEST_SUB_mybinornd( round((Counts.o2off+Counts.bg_o2off).*Counts.NBins), 0.5)./Counts.NBins;
Counts.fon_mol = TEST_SUB_mybinornd( round((Counts.o2on_mol+Counts.bg_o2on_mol).*Counts.NBins), 0.5)./Counts.NBins;
Counts.foff_mol = TEST_SUB_mybinornd( round((Counts.o2off_mol+Counts.bg_o2off_mol).*Counts.NBins), 0.5)./Counts.NBins;
Counts.fwvon = TEST_SUB_mybinornd( round((Counts.wvon+Counts.bg_wvon).*Counts.NBins), 0.5)./Counts.NBins;
Counts.fwvoff = TEST_SUB_mybinornd( round((Counts.wvoff+Counts.bg_wvoff).*Counts.NBins), 0.5)./Counts.NBins;


% subtract set f from counts to make set g
gon = (round((Counts.o2on+Counts.bg_o2on).*Counts.NBins)./Counts.NBins-Counts.fon);
goff = (round((Counts.o2off+Counts.bg_o2off).*Counts.NBins)./Counts.NBins-Counts.foff);
gon_mol = (round((Counts.o2on_mol+Counts.bg_o2on_mol).*Counts.NBins)./Counts.NBins-Counts.fon_mol);
goff_mol = (round((Counts.o2off_mol+Counts.bg_o2off_mol).*Counts.NBins)./Counts.NBins-Counts.foff_mol);
gwvon = (round((Counts.wvon+Counts.bg_wvon).*Counts.NBins)./Counts.NBins-Counts.fwvon);
gwvoff = (round((Counts.wvoff+Counts.bg_wvoff).*Counts.NBins)./Counts.NBins-Counts.fwvoff);
toc

%=====Find background of thinned profiles=====
% Counts.fon_bg = (Counts.bg_o2on.*Counts.NBins)/2;% Take mean of last data points
% Counts.fon_mol_bg = (Counts.bg_o2on_mol.*Counts.NBins)/2;% Take mean of last data points
% Counts.foff_bg = (Counts.bg_o2off.*Counts.NBins)/2;% Take mean of last data points
% Counts.foff_mol_bg = (Counts.bg_o2off_mol.*Counts.NBins)/2;% Take mean of last data points

Counts.fon_bg = (Counts.bg_o2on)/2;% Take mean of last data points
Counts.fon_mol_bg = (Counts.bg_o2on_mol)/2;% Take mean of last data points
Counts.foff_bg = (Counts.bg_o2off)/2;% Take mean of last data points
Counts.foff_mol_bg = (Counts.bg_o2off_mol)/2;% Take mean of last data points
Counts.fwvon_bg = (Counts.bg_wvon)/2;
Counts.fwvoff_bg = (Counts.bg_wvoff)/2;


% Counts.fon_bg = mean(Counts.fon(end-20:end,:),1);
% Counts.foff_bg = mean(Counts.foff(end-20:end,:),1);
% Counts.fon_mol_bg = mean(Counts.fon_mol(end-20:end,:),1);
% Counts.foff_mol_bg = mean(Counts.foff_mol(end-20:end,:),1);
% 
% Counts.gon_bg = mean(gon(end-20:end,:),1);
% Counts.goff_bg = mean(goff(end-20:end,:),1);
% Counts.gon_mol_bg = mean(gon_mol(end-20:end,:),1);
% Counts.goff_mol_bg = mean(goff_mol(end-20:end,:),1);
% 
% fwvon_bg = (Counts.bg_wvon.*Counts.NBins)/2;% Take mean of last data points
% fwvoff_bg = (Counts.bg_wvoff.*Counts.NBins)/2;% Take mean of last data points


Counts.fon = Counts.fon-Counts.fon_bg;
Counts.foff = Counts.foff-Counts.foff_bg;
Counts.fon_mol = Counts.fon_mol-Counts.fon_mol_bg;
Counts.foff_mol = Counts.foff_mol-Counts.foff_mol_bg;
Counts.fwvon = Counts.fwvon-Counts.fwvon_bg;
Counts.fwvoff = Counts.fwvoff-Counts.fwvoff_bg;

Counts.gon = gon-Counts.fon_bg;
Counts.goff = goff-Counts.foff_bg;
Counts.gon_mol = gon_mol-Counts.fon_mol_bg;
Counts.goff_mol = goff_mol-Counts.foff_mol_bg;
Counts.gwvon = gwvon-Counts.fwvon_bg;
Counts.gwvoff = gwvoff-Counts.fwvoff_bg;
% disp('creating filter')
% %=====Find optimal filters====
% [~,~,rangeWidthon,timeWidthon] = findMinE(Counts.fon,gon,fon_bg);
% 
% [~,~,rangeWidthoff,timeWidthoff] = findMinE(Counts.foff,goff,foff_bg);
% 
% [~,~,rangeWidthon_mol,timeWidthon_mol] = findMinE(Counts.fon_mol,gon_mol,fon_mol_bg);
% 
% [~,~,rangeWidthoff_mol,timeWidthoff_mol] = findMinE(Counts.foff_mol,goff_mol,foff_mol_bg);
% 
% [~,~,rangeWidthwvon,timeWidthwvon] = findMinE(Counts.fwvon,gwvon,fwvon_bg);
% [~,~,rangeWidthwvoff,timeWidthwvoff] = findMinE(Counts.fwvoff,gwvoff,fwvoff_bg);


% assign outputs to structure
% Counts.Poissonthin.timeWidthon = timeWidthon;
% Counts.Poissonthin.timeWidthoff = timeWidthoff;
% Counts.Poissonthin.timeWidthon_mol = timeWidthon_mol;
% Counts.Poissonthin.timeWidthoff_mol = timeWidthoff_mol;
% Counts.Poissonthin.rangeWidthon = rangeWidthon;
% Counts.Poissonthin.rangeWidthoff = rangeWidthoff;
% Counts.Poissonthin.rangeWidthon_mol = rangeWidthon_mol;
% Counts.Poissonthin.rangeWidthoff_mol = rangeWidthoff_mol;
% 
% Counts.Poissonthin.timeWidthwvon = timeWidthwvon;
% Counts.Poissonthin.timeWidthwvoff = timeWidthwvoff;
% Counts.Poissonthin.rangeWidthwvon = rangeWidthwvon;
% Counts.Poissonthin.rangeWidthwvoff = rangeWidthwvoff;

% Counts.Poissonthin.Ezon=Ezon;
% Counts.Poissonthin.Eton=Eton;
% Counts.Poissonthin.Ezoff=Ezoff;
% Counts.Poissonthin.Etoff=Etoff;
% Counts.Poissonthin.Ezon_mol=Ezon_mol;
% Counts.Poissonthin.Eton_mol=Eton_mol;
% Counts.Poissonthin.Ezoff_mol=Ezoff_mol;
% Counts.Poissonthin.Etoff_mol=Etoff_mol;
% 
% Counts.Poissonthin.Ezwvon=Ezwvon;
% Counts.Poissonthin.Etwvon=Etwvon;
% Counts.Poissonthin.Ezwvoff=Ezwvoff;
% Counts.Poissonthin.Etwvoff=Etwvoff;

% [Counts.o2on] = applyFilter(rangeWidthon,timeWidthon,Counts.o2on);
% [Counts.o2off] = applyFilter(rangeWidthoff,timeWidthoff,Counts.o2off);
% [Counts.o2on_mol] = applyFilter(rangeWidthon_mol,timeWidthon_mol,Counts.o2on_mol);
% [Counts.o2off_mol] = applyFilter(rangeWidthoff_mol,timeWidthoff_mol,Counts.o2off_mol);
% [Counts.wvon] = applyFilter(rangeWidthwvon,timeWidthwvon,Counts.wvon);
% [Counts.wvoff] = applyFilter(rangeWidthwvoff,timeWidthwvoff,Counts.wvoff);

end

function [counts] = applyFilter(rangeWidth,timeWidth,counts)
    %==Filter in range==
    for iii = 1:size(counts,2) %loop over time
        nz = round(4*rangeWidth(iii)); %number of grid points 
        z = (-nz:nz)';%filter grid
        kern = exp(-z.^2/rangeWidth(iii).^2); %smoothing kernel
        kern = kern/sum(sum(kern)); %normalized smoothing kernel
        counts(:,iii) = nanconv(counts(:,iii),kern,'edge','nanout'); %apply kernel to count data
    end
    %==Filter in time==
    for iii = 1:size(counts,1) %loop over range
        nz = round(4*timeWidth(iii)); %number of grid points 
        z = (-nz:nz);%filter grid
        kern = exp(-z.^2/timeWidth(iii).^2);
        kern = kern/sum(sum(kern));
        counts(:,iii) = nanconv(counts(:,iii),kern,'edge','nanout');
    end
end

function [Ez,Et,minSigz,minSigt] = findMinE(f,g,bg)
    
    %====find best filter in range====
    filt_size = logspace(-1,1.5,40); %create filter size in terms of grid points
    Ez = ones(1,length(f(1,:)),length(filt_size));
    fprintf('Range ')
    for jj = 1:length(filt_size) %loop over different filters
        fprintf('%g ',jj)
        nz = round(4*filt_size(jj)); %number of grid points 
        z = (-nz:nz)'; %filter grid
        kern = exp(-z.^2/filt_size(jj).^2); %gaussian fitler kernel in range
        if length(kern) > 1 %set gaussian filter for edge cases
            if sum(kern) == 0
                [~,it0] = min(abs(z));
                kern(it0) = 1.0;
            end
        else 
            kern = ones(1);
        end
        kern = kern/sum(sum(kern)); %normalize kernel
        fFilt = nanconv(f-bg,kern,'edge','nanout'); %apply filter
        fFilt(fFilt==0)=.001; %avoid inf in log
        Ez(:,:,jj) = sum(fFilt+bg-g.*log(fFilt+bg),1,'omitnan'); %loss function to optimize
    end
    fprintf('\n')
    [~,minEind]=min(Ez,[],3); %find minimum of loss function
    minSigz = ones(1,size(f,2));
    for ii = 1:size(f,2)
        minSigz(ii) = filt_size(minEind(ii)); %set filter width to minimum of loss function
    end

    %====find best filter in time====
    filt_size = logspace(-1.5,3,100);
    Et = ones(length(f(:,1)),1,length(filt_size));
    fprintf('time ')
    for jj = 1:length(filt_size)
        fprintf('%g ',jj)
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
    fprintf('\n')
    [~,minEind]=min(Et,[],3);
    minSigt = ones(size(f,1),1);
    for ii = 1:size(f,1)
        minSigt(ii) = filt_size(minEind(ii));
    end       
end




