clc % clears the command window
close all %removes all figures
clear
load("nysemerged.mat");
x=table2array(nyseMergedTimeTable);
[t,m]=size(x);
%% Individual stocks
stock_cum_returns=cumprod(x);
stock_set=timetable(nyseMergedTimeTable.Var1,stock_cum_returns);
plot(stock_set.Time,stock_set{:,:});
set(gca, 'YScale', 'log')
pbaspect([2 1 1])
ylabel("Wealth")
xlabel("Time")
%% bcrp performance
fn0 = @(b) (-prod(x*b));
b0=(1/m)*ones(m,1);
% No Short-selling (upper and lower bounds)
ub = ones(m,1);
lb = zeros(m,1);
% Equality constraint (fully invested)
Aeq = ones(1,m);
beq = 1;
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'Display','off');
[bcrp, fval] = fmincon(fn0,b0,[],[],Aeq,beq,lb,ub,[],options);
bcrp_cum_returns=cumprod(x*bcrp);

%% nearest neighbour K=5, L=10
load("nysemerged_nearestneighbour_results.mat");
nearestneighbour_returns=res{2}(1:end-1);
%% anticor w=20
w=20;
load("nysemerged_anticor_results.mat");
anticor_cum_returns=cumprod(anticor_returns(:,w));
%% plot
result_set=timetable(nyseMergedTimeTable.Var1,bcrp_cum_returns,'VariableNames',{'BCRP'});
result_set=[result_set timetable(nyseMergedTimeTable.Var1,nearestneighbour_returns,'VariableNames',{'Gyorfi Nearest Neighbour'})];
result_set=[result_set timetable(nyseMergedTimeTable.Var1,anticor_cum_returns,'VariableNames',{'Anticor'})];
plot(result_set.Time,result_set{:,:});
set(gca, 'YScale', 'log')
pbaspect([2 1 1])
ylabel("Wealth")
xlabel("Time")
legend(result_set.Properties.VariableNames,Location="northwest")