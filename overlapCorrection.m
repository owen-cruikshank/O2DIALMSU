%overlap correction
%%
%=====
% Overlap correction
%=====


% molecular backscatter model
[Bm,Ba,BR]= BackscatterRatiov2(ts,rm,o2off,o2off_mol,T,P,lambda_offline);

% log range corrected
mol_sig = Bm .* log(o2off_mol .* rm .^2);



log_corrected_o2on = log(o2on.*rm.^2);
log_corrected_o2off = log(o2off.*rm.^2);
log_corrected_o2off = mol_sig;

upper_O_fit = 6000;%[m]
[~,upper_O_fit_ind] = min(abs(rm-upper_O_fit));
lower_O_fit = 3000;%[m]
[~,lower_O_fit_ind] = min(abs(rm-lower_O_fit));

[~,O_fit_index] = min(abs(21-ts/60/60));

O_fit = fit(rm(lower_O_fit_ind:upper_O_fit_ind),log_corrected_o2off(lower_O_fit_ind:upper_O_fit_ind,O_fit_index),'poly1');
O_fit_slope = O_fit.p1;
O_fit_intercept = O_fit.p2;

O = log_corrected_o2off(:,O_fit_index) ./ (O_fit_intercept + rm.*O_fit_slope);

Over = O;
for i = 1:length(O)
    
    if O(i)>=1
        Over(i:end)=1;
        break;
    end
end

log_corrected_o2off_over = log_corrected_o2off ./ Over;
log_corrected_o2on_over = log_corrected_o2on ./ Over;

o2on_over = exp(log_corrected_o2on_over)./rm.^2 ;
o2off_over = exp(log_corrected_o2off_over)./rm.^2 ;

%o2on = o2on_over;
%o2off = o2off_over;

%save('fullOverlapMolchannel.mat','Over')

figure(5504)
plot(rm,o2on(:,O_fit_index));
hold on
plot(rm,o2off(:,O_fit_index));
plot(rm,o2on_over(:,O_fit_index),'--');
plot(rm,o2off_over(:,O_fit_index),'--');
%plot(rm,O_fit_intercept + rm.*O_fit_slope);
hold off
title('log range corrected counts')
legend('on','off','on over','off over')

figure(5505)
plot(rm,log_corrected_o2on(:,O_fit_index));
hold on
plot(rm,log_corrected_o2off(:,O_fit_index));
plot(rm,log_corrected_o2off_over(:,O_fit_index),'--');
plot(rm,log_corrected_o2on_over(:,O_fit_index),'--');
plot(rm,O_fit_intercept + rm.*O_fit_slope,'.-');
hold off
title('log range corrected counts')

figure(5506)
plot(rm,O);
hold on
plot(rm,rm./rm)
plot(rm,Over)
title('overlap function')
hold off