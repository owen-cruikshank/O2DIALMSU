function [ res ] = TEST_SUB_mybinornd( N, p )
    [row_cnt, col_cnt] = size(N);
    res = zeros(row_cnt, col_cnt);
    for ii=1:row_cnt
       for jj=1:col_cnt
           if isnan(N(ii,jj))
               res(ii,jj)=NaN;
           else
               res(ii, jj) = sum(rand(1,N(ii,jj))<p);

               %res(ii, jj) = sum(randi([0 1],N(ii,jj),1));
           end
       end
    end
end