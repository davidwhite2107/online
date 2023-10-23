clc % clears the command window
close all %removes all figures
clear
%% 0. Synthetic data
% 7 stocks (M) and 500 ticks (n) normal distributed with mean drift 1%
m = 36;
n = 500;
k0 = 21;
% 7 random stock returns
rng(1,'twister');
x = 0.03*randn(n,m)+0.01;
%prices
p1 = ret2price(x);
p1=p1(2:end,:)./p1(1:end-1,:);
t = m;
r = tick2ret(p1,[],'Simple');
%% anticor
w=15;
b_mine=anticor(p1,w);
b_tim=anticor_tim(p1,w,true);