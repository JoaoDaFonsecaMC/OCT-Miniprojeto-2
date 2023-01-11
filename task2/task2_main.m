%MP2 - task II
clear all;
format long;

%parameters
N = 100;
t_0 = 0;
t_f = 1;

%time step
h = (t_f-t_0)/N;

time = t_0:h:t_f;

X0 = ones(N+1,5);

options=optimoptions('fmincon','MaxFunctionEvaluations',10000000000000,'MaxIterations',10000000000000);
[X,fval,exitflag,output] = fmincon(@cost_function,X0,[],[],[],[],[],[],@restrictions,options);

plot(time,X(:,1));