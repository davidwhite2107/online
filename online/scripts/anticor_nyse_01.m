clc % clears the command window
close all %removes all figures
clear

%% anticor on NYSE dataset
load("nyse.mat");
x=table2array(nyseTimeTable);
[t0,m]=size(x);
w=20;
b=anticor(x, w);
anticor_returns=[];
for t=1:length(b)
  anticor_returns(t,:) = (b(t,:) * x(t,:)');
end
ANTICOR=cumprod(anticor_returns);
port_ret=ANTICOR(end);