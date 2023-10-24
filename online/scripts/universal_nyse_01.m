%% universal algorithm on pairs of stocks from NYSE dataset
clc % clears the command window
close all %removes all figures
clear
load("nyse.mat");

%iroqu_kinar
iroqu_kinar_tickers={'iroqu', 'kinar'};
iroqu_kinar=nyseTimeTable{:,iroqu_kinar_tickers};
[t,m]=size(iroqu_kinar);

%UP
iroqu_kinar_b=universal(iroqu_kinar);
iroqu_kinar_returns=[];
%accumulated returns
for t=1:length(iroqu_kinar_b)
  iroqu_kinar_returns(t,:) = (iroqu_kinar_b(t,:) * iroqu_kinar(t,:)');
end
IK_return=cumprod(iroqu_kinar_returns);

%individual stocks
iroqu=cumprod(nyseTimeTable{:,{'iroqu'}});
kinar=cumprod(nyseTimeTable{:,{'kinar'}});


%bcrp performance
fn0 = @(b) (-prod(iroqu_kinar*b));
b0=(1/m)*ones(m,1);
% No Short-selling (upper and lower bounds)
ub = ones(m,1);
lb = zeros(m,1);
% Equality constraint (fully invested)
Aeq = ones(1,m);
beq = 1;
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'Display','off');
[IK_bcrp, fval] = fmincon(fn0,b0,[],[],Aeq,beq,lb,ub,[],options);
IK_bcrp_returns=cumprod(iroqu_kinar*IK_bcrp);


IK_set=timetable(nyseTimeTable.Var1,IK_return,'VariableNames',{'Universal Portfolio'});
IK_set=[IK_set timetable(nyseTimeTable.Var1,iroqu,'VariableNames',{'IROQU'})];
IK_set=[IK_set timetable(nyseTimeTable.Var1,kinar,'VariableNames',{'KINAR'})];
IK_set=[IK_set timetable(nyseTimeTable.Var1,IK_bcrp_returns,'VariableNames',{'BCRP'})];
plot(IK_set.Time,IK_set{:,:});
ylabel("Wealth")
xlabel("Time")
title("IROQU & KINAR (NYSE)")
legend(IK_set.Properties.VariableNames,Location="northwest")
%%
%comme_meico
comme_meico_tickers={'comme', 'meico'};
comme_meico=nyseTimeTable{:,comme_meico_tickers};
[t,m]=size(comme_meico);

%UP
comme_meico_b=universal(comme_meico);
comme_meico_returns=[];
for t=1:length(comme_meico_b)
  comme_meico_returns(t,:) = (comme_meico_b(t,:) * comme_meico(t,:)');
end
CM_return=cumprod(comme_meico_returns);

%individual stocks
comme=cumprod(nyseTimeTable{:,{'comme'}});
meico=cumprod(nyseTimeTable{:,{'meico'}});

%bcrp performance
fn1 = @(b) (-prod(comme_meico*b));
[CM_bcrp, fval] = fmincon(fn1,b0,[],[],Aeq,beq,lb,ub,[],options);
CM_bcrp_returns=cumprod(comme_meico*CM_bcrp);

%plot returns
CM_set=timetable(nyseTimeTable.Var1,CM_return,'VariableNames',{'Universal Portfolio'});
CM_set=[CM_set timetable(nyseTimeTable.Var1,comme,'VariableNames',{'COMME'})];
CM_set=[CM_set timetable(nyseTimeTable.Var1,meico,'VariableNames',{'MEICO'})];
CM_set=[CM_set timetable(nyseTimeTable.Var1,CM_bcrp_returns,'VariableNames',{'BCRP'})];
plot(CM_set.Time,CM_set{:,:});
ylabel("Wealth")
xlabel("Time")
title("COMME & MEICO (NYSE)")
legend(CM_set.Properties.VariableNames,Location="northwest")
%%
%comme_kinar
comme_kinar_tickers={'comme', 'kinar'};
comme_kinar=nyseTimeTable{:,comme_kinar_tickers};
[t,m]=size(comme_kinar);

%UP
comme_kinar_b=universal(comme_kinar);
comme_kinar_returns=[];
for t=1:length(comme_kinar_b)
  comme_kinar_returns(t,:) = (comme_kinar_b(t,:) * comme_kinar(t,:)');
end
CK_return=cumprod(comme_kinar_returns);

%individual stocks
comme=cumprod(nyseTimeTable{:,{'comme'}});
meico=cumprod(nyseTimeTable{:,{'kinar'}});

%bcrp performance
fn2 = @(b) (-prod(comme_kinar*b));
[CK_bcrp, fval] = fmincon(fn2,b0,[],[],Aeq,beq,lb,ub,[],options);
CK_bcrp_returns=cumprod(comme_kinar*CK_bcrp);

%plot returns
CK_set=timetable(nyseTimeTable.Var1,CK_return,'VariableNames',{'Universal Portfolio'});
CK_set=[CK_set timetable(nyseTimeTable.Var1,comme,'VariableNames',{'COMME'})];
CK_set=[CK_set timetable(nyseTimeTable.Var1,kinar,'VariableNames',{'KINAR'})];
CK_set=[CK_set timetable(nyseTimeTable.Var1,CK_bcrp_returns,'VariableNames',{'BCRP'})];
plot(CK_set.Time,CK_set{:,:});
ylabel("Wealth")
xlabel("Time")
title("COMME & KINAR (NYSE)")
legend(CK_set.Properties.VariableNames,Location="northwest")
%%
%ibm_coke
ibm_coke_tickers={'ibm', 'coke'};
ibm_coke=nyseTimeTable{:,ibm_coke_tickers};
[t,m]=size(ibm_coke);

%UP
ibm_coke_b=universal(ibm_coke);
ibm_coke_returns=[];
for t=1:length(ibm_coke_b)
  ibm_coke_returns(t,:) = (ibm_coke_b(t,:) * ibm_coke(t,:)');
end
IC_return=cumprod(ibm_coke_returns);

%individual stocks
ibm=cumprod(nyseTimeTable{:,{'ibm'}});
coke=cumprod(nyseTimeTable{:,{'coke'}});

%bcrp performance
fn3 = @(b) (-prod(ibm_coke*b));
[IC_bcrp, fval] = fmincon(fn3,b0,[],[],Aeq,beq,lb,ub,[],options);
IC_bcrp_returns=cumprod(comme_kinar*IC_bcrp);

%plot returns
IC_set=timetable(nyseTimeTable.Var1,IC_return,'VariableNames',{'Universal Portfolio'});
IC_set=[IC_set timetable(nyseTimeTable.Var1,ibm,'VariableNames',{'IBM'})];
IC_set=[IC_set timetable(nyseTimeTable.Var1,coke,'VariableNames',{'COKE'})];
IC_set=[IC_set timetable(nyseTimeTable.Var1,IC_bcrp_returns,'VariableNames',{'BCRP'})];
plot(IC_set.Time,IC_set{:,:});
ylabel("Wealth")
xlabel("Time")
title("IBM & COKE (NYSE)")
legend(IC_set.Properties.VariableNames,Location="northwest")
