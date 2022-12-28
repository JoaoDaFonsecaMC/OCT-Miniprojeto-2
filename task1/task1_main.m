%MP 2, Task I
clear all;
clc;

%minimization options
options = optimset('LargeScale','off','Display','iter');

%- Problem 1 -

%initial guess
x0 = [0 0 0 0];
%result
[x,fval1,exitflag1,output1] = fmincon(@objfun1,x0,[],[],[],[],[],[],@confun1,options);

%- Problem 2 -

%initial guess
y0= [0 0 0 0 0 0 0];
%result
[y,fval2,exitflag2,output2] = fmincon(@objfun2,y0,[],[],[],[],[],[],@confun2,options);