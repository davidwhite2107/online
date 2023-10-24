%% nearest neighbour algorithm on NYSE dataset
clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
x=table2array(nyseTimeTable);
K=5;
L=10;
%%
[b, S, h, SH]=nearestneighbour(x, 1:K, 1:L);
nn={b, S, h, SH};
save("nyse_nearestneighbour_results", "nn");
