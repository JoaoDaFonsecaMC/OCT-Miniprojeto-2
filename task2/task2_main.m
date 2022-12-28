%MP2 - task II
clear all;
clc;

phi0 = [0];
options = optimset('LargeScale','off','Display','iter');
[phi,fval,exitflag,output] = fmincon(@cost_function,phi0,[],[],[],[],-0.5,0.5,[],options);