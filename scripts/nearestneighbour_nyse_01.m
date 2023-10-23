clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
x=table2array(nyseTimeTable);
K=5;
L=10;
%%
[b, S, h, SH]=nearestneighbour(x, 1:K, 1:L);
%%
[t0, m0]=size(x);
%%
nn_returns=zeros(t0, 1);
for t=1:t0
  nn_returns(t,:) = (b(t,:) * x(t,:)');
end
NN=cumprod(nn_returns(1:t0,:));
port_ret=NN(end);
%%
x2=reshape(x',m0,1,t0);
p=gyorfinn(x2,1:5,1:10,[],'gyorfi_opt');
p=offline(p);