clc % clears the command window
close all %removes all figures
clear
load("nyse_anticor_results");
anticor_returns=anticor_returns(:,2:end);
%% psr and dsr
%returns of the 'strategies'
returns=anticor_returns-1;
%num strategies
[t,m]=size(returns);
n=m;
%euler-mascheroni constant
gamma=0.5772;
%1. compute sharpe ratios, skewness and kurtosis for each strategy
%mean and variance of each strategy's returns
mu=mean(returns, 1);
sigma=var(returns, 1);
%skewness and kurtosis of each strategy's returns
skew=skewness(returns, 1);
kurt=kurtosis(returns, 1);
%sharpe ratios for each strategy's returns
sr=(mu./sqrt(sigma))';

%2. find optimal clusters of strategies based on sharpe ratio
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton','replicate',5));
eva = evalclusters(sr,myfunc,'CalinskiHarabasz','klist',1:10);
k=eva.OptimalK;
% make the clustering reproducible
rng(1, "twister");
clusters=kmeans(sr, k);

%3. find average sharpe ratio, skewness and kurtosis in each cluster
ave_strategy_srs=zeros(k,1);
ave_strategy_skew=zeros(k,1);
ave_strategy_kurt=zeros(k,1);
for i=1:k
    ave_strategy_srs(i)=mean(sr(clusters==i));
    ave_strategy_skew(i)=mean(skew(clusters==i));
    ave_strategy_kurt(i)=mean(kurt(clusters==i));
end
%4. compute benchmark sharpe ratio
var_sr_k=var(ave_strategy_srs);
sr_0=sqrt(var_sr_k)*((1-gamma)*norminv(1-(1/k))+gamma*norminv(1-(1/k)*(1/exp(1))));
%5. use benchmark sharpe ratio and compute dsr and psr
psr=zeros(k,1);
dsr=zeros(k, 1);
for i=1:k
    %calculate psr for strategy i
    psr(i) = ((ave_strategy_srs(i)-sr_0)*sqrt(k-1))/(sqrt(1-ave_strategy_srs(i)*ave_strategy_skew(i)+ave_strategy_srs(i)^2*((ave_strategy_kurt(i)-1)/4)));
    %calculate dsr for strategy i
    dsr(i)=normcdf(psr(i));
end
%% pbo
%1. partition returns matrix, M, into M_s matrices
%pick size of partition
s=16;
index=1:s;
%size of each m matrix
size_m=floor(t/s);
%cell array to store the matrices
m_s={};
%partition into m_s
for i=0:s-1
    start=(i*size_m)+1;
    finish=(i+1)*size_m;
    if i==s-1
        finish=t;
    end
    m_s{i+1}=returns(start:finish,:);
end
%2. Form all combinations of M_s, taken in groups of S/2, call each
%combination c
C_indices=nchoosek(index, s/2);
length_C=height(C_indices);
lambda_c=zeros(length_C, 1);
%3. For each c,
for i=1:length_C
    %3.a. Form training set J by taking S/2 M_s of c, preserving their order
    %find index for combination c
    j_index=C_indices(i,:);
    %complement of c in index to get Jbar index
    jbar_index=index(setdiff(1:end,j_index));
    J=vertcat(m_s{j_index});
    %3.b. Form the test set Jbar, by taking complement of J in M
    Jbar=vertcat(m_s{jbar_index});
    %3.c.i.  Form vector of Sharpe ratios R^c for each column in J
    mu=mean(J, 1);
    sigma=var(J, 1);
    sr_J=mu./sqrt(sigma);
    %3.c.ii. Rank Sharpe ratios in R^c and store as r_c
    %3.d. Repeat 3.c. for J_bar
    %3.d.i.  Form vector of Sharpe ratios R^c for each column in Jbar
    mu_bar=mean(Jbar, 1);
    sigma_bar=var(Jbar, 1);
    sr_Jbar=mu_bar./sqrt(sigma_bar);
    %3.d.ii. Rank Sharpe ratios in R^c and store as r_c
    %3.e. Find the best performing strategy in-sample, r^c_n*, and its index, n*
    [max_sr_J, nstar]=max(sr_J);
    %3.f. Find relative rank out-of-sample, wbar^c = rbar^c_n*/(n+1)
    sr_nstar=sr_Jbar(:,nstar);
    rbar_nstar=find(sort(sr_Jbar)==sr_nstar);
    wbar=rbar_nstar/(n+1);
    %3.g. Define lambda=logit(wbar^c)
    lambda_c(i,:)=log((wbar)/(1-wbar));
end
%4. Find relative frequency of each lambda to create a distribution
%f(lambda)
f_lambda=lambda_c(lambda_c<=0,:);
%5. Find PBO: phi=int_{-Inf}^0 f(lambda) dlambda
pbo=length(f_lambda)/length_C;
%%
histogram(lambda_c, 50, 'Normalization', 'probability')
xlabel("Logits")
ylabel("Probability")
xline(0, '--', 'color', 'red')