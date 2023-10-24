%% performance of all algorithms on iroqu kinar from NYSE dataset
clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");
iroqu_kinar_tickers={'iroqu', 'kinar'};
x=nyseTimeTable{:,iroqu_kinar_tickers};
[t,m]=size(x);

%% individual stocks
iroqu=cumprod(nyseTimeTable{:, {'iroqu'}});
kinar=cumprod(nyseTimeTable{:, {'kinar'}});
stock_set=timetable(nyseTimeTable.Var1,iroqu,'VariableNames',{'iroqu'});
stock_set=[stock_set timetable(nyseTimeTable.Var1,kinar,'VariableNames',{'kinar'})];
plot(stock_set.Time,stock_set{:,:});
pbaspect([2 1 1])
ylabel("Wealth")
xlabel("Time")
legend(stock_set.Properties.VariableNames,Location="northeast")
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
load("iroqu_kinar_nn");
iroqu_kinar_S=iroqu_kinar_S(1:end-1);
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
result_set=timetable(nyseTimeTable.Var1,bcrp_returns,'VariableNames',{'BCRP'});
result_set=[result_set timetable(nyseTimeTable.Var1,iroqu_kinar_S,'VariableNames',{'Gyorfi Nearest Neighbour'})];
result_set=[result_set timetable(nyseTimeTable.Var1,anticor_return,'VariableNames',{'Anticor'})];
result_set=[result_set timetable(nyseTimeTable.Var1,up_return,'VariableNames',{'Universal Portfolio'})];
plot(result_set.Time,result_set{:,:});
set(gca, 'YScale', 'log')
pbaspect([2 1 1])
ylabel("Wealth")
xlabel("Time")
legend(result_set.Properties.VariableNames,Location="northwest")