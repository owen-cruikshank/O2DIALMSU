function [data,bins] = intSum(inData,inDatats,NBins,ts)
%Integrate data to new time grid

%ins
%inData = [rangextime] input data
%inDatats = [1xtime] input data time grid
%NBins = [1xtime] number of integrated bins already

%outs
%data = [rangextime] integrted and averaed data
%bins = [1xtime] number of averaged bins

increment = 1;
bins = zeros(1,size(ts,2));
binsNew = zeros(1,size(ts,2));
data = zeros(size(inData,1),size(ts,2));
for ii = 1:size(inDatats,2)
    if increment > size(ts,2)%Check if increment is too long
        data(:,increment-1) = data(:,increment-1)+inData(:,ii);
        bins(:,increment-1) = bins(:,increment-1)+NBins(:,ii);
        binsNew(:,increment-1) = binsNew(:,increment-1)+1;
        
    elseif inDatats(ii)<=ts(increment)%check if you can continue to add
        data(:,increment) = data(:,increment)+inData(:,ii);
        bins(:,increment) = bins(:,increment)+NBins(:,ii);
        binsNew(:,increment) = binsNew(:,increment)+1;
    else
        increment = increment+1;
        if increment > size(ts,2)%check if done
            
        else %add onto new bin
            data(:,increment) = data(:,increment)+inData(:,ii);
            bins(:,increment) = bins(:,increment)+NBins(:,ii);
            binsNew(:,increment) = binsNew(:,increment)+1;
        end

    end
end
%average
data = data./binsNew;
        