clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
x=table2array(nyseTimeTable);
[t0,m0]=size(x);
%%
w=70;
anticor_returns=zeros(t0,w);
anticor_cum_ret=zeros(t0,w);
for w0=2:w
    b=anticor(x, w0);
    for t=1:t0
        anticor_returns(t,w0) = (b(t,:) * x(t,:)');
    end
    anticor_cum_ret(:,w0)=cumprod(anticor_returns(:,w0));
end
window_returns=anticor_cum_ret(t0,2:w);
save('nyse_anticor_results', 'anticor_returns');
%%
w=70;
load('nyse_anticor_results.mat');
anticor_cum_ret=cumprod(anticor_returns);
window_returns=anticor_cum_ret(t0,2:w);
[~,ind]=max(window_returns);
w_opt=ind+1;
best_stock=cumprod(x(:,30));
best_stock=best_stock(end);
windows=2:w;
plot(windows,window_returns(end,:), 'black', 'LineWidth', 1.5);
set(gca, 'YScale', 'log')
ylabel("Wealth")
xlabel("Window Size (w)")
xlim([2 w])
xticks( [2, xticks()] )
xline(w_opt, '--', 'color', 'blue')
yline(best_stock, '--', 'color','r')
title("NYSE: Anticor_w vs. window size")
legend('Anticor_w', 'Optimal w', 'Best Stock',Location="northeast")