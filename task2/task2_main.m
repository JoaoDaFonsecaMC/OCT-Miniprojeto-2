%MP2 - task II
clear all;
format long;

%parameters
N = 100;

X0 = ones(5,N+1);

options=optimoptions('fmincon','MaxFunctionEvaluations',10000000000000,'MaxIterations',10000000000000);
[X,fval,exitflag,output] = fmincon(@cost_function,X0,[],[],[],[],[],[],@restrictions,options);

plot(time,X(1,:));