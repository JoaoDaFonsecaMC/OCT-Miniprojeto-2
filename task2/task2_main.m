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

%initial guess
X0 = ones(N+1,5);

%minimization
options=optimoptions('fmincon','MaxFunctionEvaluations',10000000000000,'MaxIterations',10000000000000);
[X,fval,exitflag,output] = fmincon(@cost_function,X0,[],[],[],[],[],[],@restrictions,options);

%plots
figure()
plot(time,X(:,1));
title("Phi vs time");
xlabel ('time[s]');
ylabel ('phi[rad]');

figure()
plot(time,X(:,2));
title("u vs time");
xlabel ('time[s]');
ylabel ('u[km/s]');

figure()
plot(time,X(:,3));
title("v vs time");
xlabel ('time[s]');
ylabel ('v[km/s]');

figure()
plot(time,X(:,4));
title("r vs time");
xlabel ('time[s]');
ylabel ('r[km]');

figure()
plot(time,X(:,5));
title("theta vs time");
xlabel ('time[s]');
ylabel ('theta[rad]');