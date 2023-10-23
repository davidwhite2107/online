clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
x=table2array(nyseTimeTable);
[t0, m0]=size(x);
%%
K1=10;
L1=20;
[b1, S1, h1, SH1]=nearestneighbour(x, 1:K1, 1:L1);
K2=15;
L2=30;
[b2, S2, h2, SH2]=nearestneighbour(x, 1:K2, 1:L2);