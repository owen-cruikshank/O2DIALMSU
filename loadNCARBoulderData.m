function [Data] = loadNCARBoulderData(spanDays,path)

% for ii = 1:length(spanDays)
ii=1;
    fullDataPath = strcat(path,'mpd05.',num2str(year(spanDays(ii))),sprintf('%02d',month(spanDays(ii))),sprintf('%02d',day(spanDays(ii))),'.MatlabPreload.mat');
    load(fullDataPath)
    
    dataNames = fieldnames(Data.TimeSeries);
    for jj = 1:length(dataNames)
        Data.(dataNames{jj})=Data.TimeSeries.(dataNames{jj});
    end
    
%     if ii==1
%         DataOut = Data;
%     end
    

end

% function [] = recursiveConcatination(DataOut,Data,numb)
%     
%     names = fieldnames(Data);
%     if numb = length(
%     
%     elseif ~isstruct(Data.(names{1}))
%         for jj = 1:length(names)
%             if strcmp('TimeStamp',names{jj})
%                 DataOut.(names{1}) = [DataOut.(names{jj}) Data.(names{jj})+DataOut.(names{jj})];
%             else
%                 DataOut.(names{1}) = [DataOut.(names{jj}) Data.(names{jj})];
%             end
%         end
%     else
%         recursiveConcatination(DataOut(names{1}),Data(names{1},numb)
%     end
% end

