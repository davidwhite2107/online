clc % clears the command window
close all %removes all figures
clear
load("b_fayyaaz.mat");
%% 0. Synthetic data
% 7 stocks (M) and 500 ticks (n) normal distributed with mean drift 1%
m = 7;
n = 500;
k0 = 21;
% 7 random stock returns
rng(1,'twister');
x = 0.03*randn(n,m)+0.01;
%prices
p1 = ret2price(x);
p1=p1(2:end,:)./p1(1:end-1,:);
K = 5;
L=10;
t = m+1;
r = tick2ret(p1,[],'Simple');
p2=reshape(p1',m,1,n);
%% function
[b, S, h, SH]=nearestneighbour(p1, 1:K, 1:L);
%%
ret=zeros(n,1);
for t=1:n
    ret(t)=b(t,:)*p1(t,:)';
end
port=cumprod(ret);