clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
x=table2array(nyseTimeTable);
[t0,m0]=size(x);
K=5;
L=10;

%% iroqu_kinar
iroqu_kinar_tickers={'iroqu', 'kinar'};
iroqu_kinar=nyseTimeTable{:,iroqu_kinar_tickers};
[iroqu_kinar_b, iroqu_kinar_S, iroqu_kinar_h, iroqu_kinar_SH]=nearestneighbour(iroqu_kinar, 1:K, 1:L);
save("iroqu_kinar_nn", "iroqu_kinar_S");

%% comme_meico
comme_meico_tickers={'comme', 'meico'};
comme_meico=nyseTimeTable{:,comme_meico_tickers};
[comme_meico_b, comme_meico_S, comme_meico_h, comme_meico_SH]=nearestneighbour(comme_meico, 1:K, 1:L);
comme_meico={comme_meico_b, comme_meico_S, comme_meico_h, comme_meico_SH};
save("comme_meico_nn", "comme_meico");

%% comme_kinar
comme_kinar_tickers={'comme', 'kinar'};
comme_kinar=nyseTimeTable{:,comme_kinar_tickers};
[comme_kinar_b, comme_kinar_S, comme_kinar_h, comme_kinar_SH]=nearestneighbour(comme_kinar, 1:K, 1:L);
comme_kinar={comme_kinar_b, comme_kinar_S, comme_kinar_h, comme_kinar_SH};
save("comme_kinar_nn", "comme_kinar");

%% ibm_coke
ibm_coke_tickers={'ibm', 'coke'};
ibm_coke=nyseTimeTable{:,ibm_coke_tickers};
[ibm_coke_b, ibm_coke_S, ibm_coke_h, ibm_coke_SH]=nearestneighbour(ibm_coke, 1:K, 1:L);
ibm_coke={ibm_coke_b, ibm_coke_S, ibm_coke_h, ibm_coke_SH};
save("nyse_ibm_coke_nn", "ibm_coke");

%%
load("nysemerged.mat");
merged_ibm_coke=nyseMergedTimeTable{:,ibm_coke_tickers};
[merged_ibm_coke_b, merged_ibm_coke_S, merged_ibm_coke_h, merged_ibm_coke_SH]=nearestneighbour(merged_ibm_coke, 1:K, 1:L);
nysemerged_ibm_coke={merged_ibm_coke_b, merged_ibm_coke_S, merged_ibm_coke_h, merged_ibm_coke_SH};
save("nysemerged_ibm_coke_nn", "nysemerged_ibm_coke");