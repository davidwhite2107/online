%% anticor on pairs of stocks from the NYSE dataset
clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
x=table2array(nyseTimeTable);
[t0,m0]=size(x);
%%
w=20;

%iroqu_kinar
iroqu_kinar_tickers={'iroqu', 'kinar'};
iroqu_kinar=nyseTimeTable{:,iroqu_kinar_tickers};
[t,m]=size(iroqu_kinar);

%anticor
iroqu_kinar_b=anticor(iroqu_kinar, w);
iroqu_kinar_returns=zeros(t,1);
%accumulated returns
for t=1:length(iroqu_kinar_b)
  iroqu_kinar_returns(t,:) = (iroqu_kinar_b(t,:) * iroqu_kinar(t,:)');
end
IK_return=cumprod(iroqu_kinar_returns);

%%
%comme_meico
comme_meico_tickers={'comme', 'meico'};
comme_meico=nyseTimeTable{:,comme_meico_tickers};
[t,m]=size(comme_meico);

%UP
comme_meico_b=anticor(comme_meico, w);
comme_meico_returns=zeros(t, 1);
for t=1:length(comme_meico_b)
  comme_meico_returns(t,:) = (comme_meico_b(t,:) * comme_meico(t,:)');
end
CM_return=cumprod(comme_meico_returns);

%%
%comme_kinar
comme_kinar_tickers={'comme', 'kinar'};
comme_kinar=nyseTimeTable{:,comme_kinar_tickers};
[t,m]=size(comme_kinar);

%UP
comme_kinar_b=anticor(comme_kinar, w);
comme_kinar_returns=zeros(t, 1);
for t=1:length(comme_kinar_b)
  comme_kinar_returns(t,:) = (comme_kinar_b(t,:) * comme_kinar(t,:)');
end
CK_return=cumprod(comme_kinar_returns);
%%
%ibm_coke
ibm_coke_tickers={'ibm', 'coke'};
ibm_coke=nyseTimeTable{:,ibm_coke_tickers};
[t,m]=size(ibm_coke);

%UP
ibm_coke_b=anticor(ibm_coke, w);
ibm_coke_returns=zeros(t, 1);
for t=1:length(ibm_coke_b)
  ibm_coke_returns(t,:) = (ibm_coke_b(t,:) * ibm_coke(t,:)');
end
IC_return=cumprod(ibm_coke_returns);