function [ res ] = TEST_SUB_mybinornd2( N, p )
    %Custom fast binomial 
    [row_cnt, col_cnt] = size(N);
    res = zeros(row_cnt, col_cnt);
   % for ii=1:row_cnt
% %        for jj=1:col_cnt
% %            if isnan(N(1,jj))
% %                res(:,jj)=NaN;
% %            else
% %                maxN = max(N(:,jj));
% %                res(:, jj) = sum(rand(row_cnt,maxN)<p,2);
% %            end
% %       % end
% %        end

              for jj=1:col_cnt
           if isnan(N(1,jj))
               res(:,jj)=NaN;
           else
               maxN = max(N(:,jj));
               res(:, jj) = sum(rand(row_cnt,maxN)<p,2);
           end
      % end
    end

% for jj=1:col_cnt
%     %for ii=1:row_cnt
%        
%            if isnan(N(1,jj))
%                res(:,jj)=NaN;
%            else
%                aaa = max(N(:,jj));
%                res(:, jj) = sum(rand(row_cnt,aaa)<p,2);
%            end
%       % end
% end


end