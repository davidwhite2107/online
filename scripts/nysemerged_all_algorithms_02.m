clc % clears the command window
close all %removes all figures
clear
load("nysemerged.mat");
ibm_coke_tickers={'ibm', 'coke'};
x=nyseMergedTimeTable{:,ibm_coke_tickers};
[t,m]=size(x);

%% individual stocks
ibm=cumprod(nyseMergedTimeTable{:, {'ibm'}});
coke=cumprod(nyseMergedTimeTable{:, {'coke'}});
stock_set=timetable(nyseMergedTimeTable.Var1,ibm,'VariableNames',{'ibm'});
stock_set=[stock_set timetable(nyseMergedTimeTable.Var1,coke,'VariableNames',{'coke'})];
plot(stock_set.Time,stock_set{:,:});
pbaspect([2 1 1])
ylabel("Wealth")
xlabel("Time")
legend(stock_set.Properties.VariableNames,Location="northwest")
%% anticor
w=20;
anticor_b=anticor(x, w);
anticor_returns=zeros(t, 1);
%accumulated returns
for t=1:length(anticor_b)
  anticor_returns(t,:) = (anticor_b(t,:) * x(t,:)');
end
anticor_return=cumprod(anticor_returns);

%% universal
up_b=universal(x);
up_returns=zeros(t, 1);
%accumulated returns
for t=1:length(up_b)
  up_returns(t,:) = (up_b(t,:) * x(t,:)');
end
up_return=cumprod(up_returns);

%% nearest neighbour
load("nysemerged_ibm_coke_nn");
%% bcrp
%bcrp performance
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
bcrp_returns=cumprod(x*bcrp);
%%
result_set=timetable(nyseMergedTimeTable.Var1,bcrp_returns,'VariableNames',{'BCRP'});
%result_set=[result_set timetable(nyseMergedTimeTable.Var1,iroqu_kinar_S,'VariableNames',{'Gyorfi Nearest Neighbour'})];
result_set=[result_set timetable(nyseMergedTimeTable.Var1,anticor_return,'VariableNames',{'Anticor'})];
result_set=[result_set timetable(nyseMergedTimeTable.Var1,up_return,'VariableNames',{'Universal Portfolio'})];
plot(result_set.Time,result_set{:,:});
set(gca, 'YScale', 'log')
pbaspect([2 1 1])
ylabel("Wealth")
xlabel("Time")
legend(result_set.Properties.VariableNames,Location="northwest")