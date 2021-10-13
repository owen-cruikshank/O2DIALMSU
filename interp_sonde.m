%  Putting sonda data on the same range spacing for use with DIAL
%  Temperature Performance modeling programs

function [T_int,P_int,rm_int] = interp_sonde(T, P, r_TP, del_r)

%  ============= inputs
rangebin = del_r;


alt_0 = r_TP(1);
alt = r_TP;
range = alt-alt_0;
i_range = length(range);

tk = T;
pre = P;


j = 1;
for i = 1:i_range
    r = range(i);
    if r > j*rangebin
       p = i; 
       dr1 = range(i)-range(i-1);
       dr2 = j*rangebin - range(i-1);
       
       dt = tk(i) - tk(i-1);
       T_int(j) = tk(i-1) + dt*dr2/dr1;
       
       dp = pre(i) - pre(i-1);
       P_int(j) = pre(i-1) + dp*dr2/dr1;
       
       
       rm_int(j) = j*rangebin;
       
       j=j+1;
    end
end

