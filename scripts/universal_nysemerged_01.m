clc % clears the command window
close all %removes all figures
clear
load("nysemerged.mat");

%% ibm_coke
ibm_coke_tickers={'ibm', 'coke'};
ibm_coke=nyseMergedTimeTable{:,ibm_coke_tickers};
[t,m]=size(ibm_coke);

%UP
ibm_coke_b=universal(ibm_coke);
ibm_coke_returns=[];
for t=1:length(ibm_coke_b)
  ibm_coke_returns(t,:) = (ibm_coke_b(t,:) * ibm_coke(t,:)');
end
IC_return=cumprod(ibm_coke_returns);

%individual stocks
ibm=cumprod(nyseMergedTimeTable{:,{'ibm'}});
coke=cumprod(nyseMergedTimeTable{:,{'coke'}});

%bcrp performance
fn3 = @(b) (-prod(ibm_coke*b));
b0=(1/m)*ones(m,1);
% No Short-selling (upper and lower bounds)
ub = ones(m,1);
lb = zeros(m,1);
% Equality constraint (fully invested)
Aeq = ones(1,m);
beq = 1;
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'Display','off');
[IC_bcrp, fval] = fmincon(fn3,b0,[],[],Aeq,beq,lb,ub,[],options);
IC_bcrp_returns=cumprod(ibm_coke*IC_bcrp);

%plot returns
IC_set=timetable(nyseMergedTimeTable.Var1,IC_return,'VariableNames',{'UP'});
IC_set=[IC_set timetable(nyseMergedTimeTable.Var1,ibm,'VariableNames',{'IBM'})];
IC_set=[IC_set timetable(nyseMergedTimeTable.Var1,coke,'VariableNames',{'COKE'})];
IC_set=[IC_set timetable(nyseMergedTimeTable.Var1,IC_bcrp_returns,'VariableNames',{'BCRP'})];
plot(IC_set.Time,IC_set{:,:});
ylabel("Wealth")
xlabel("Time")
title("IBM & COKE (NYSE)")
legend(IC_set.Properties.VariableNames,Location="northwest")

%%