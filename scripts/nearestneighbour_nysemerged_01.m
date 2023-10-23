clc % clears the command window
close all %removes all figures
clear
load("nysemerged.mat");
x=table2array(nyseMergedTimeTable);
K=5;
L=10;
[b, S, h, SH]=nearestneighbour(x, 1:K, 1:L);
res={b, S, h, SH};
save("nysemerged_nearestneighbour_results.mat", "res");