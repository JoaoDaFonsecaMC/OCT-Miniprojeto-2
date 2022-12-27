%MP 3, Task I
clear all;
clc;

x0 = [-1 1];
options = optimset('LargeScale','off','Display','iter');
[x,fval,exitflag,output] = fmincon(@objfun,x0,[],[],[],[],[],[],@confun,options);

% Objective Function (M-file)
function f = objfun(x)
% objective function
f= x(1)^(2)+x(2)^(2)+2*x(3)^(2)+x(4)^(2)-5*x(1)-5*x(2)-21*x(3)+7*x(4);
end
% Constraints Function (M-file)
function [c, ceq] = confun(x)
% Nonlinear inequality constraints:
c = [-8+x(1)-x(2)+x(3)-x(4)+x(1)^(2)+x(2)^(2)+x(3)^(2)+x(4)^(2); 
 -10-x(1)-x(4)+x(1)^(2)+2*x(2)^(2)+x(3)^(2)+2*x(4)^(2)];
% no nonlinear equality constraints:
ceq = [5-2*x(1)+x(2)+x(4)-2*x(1)^(2)-x(2)^(2)-2*x(3)^(2)];
end