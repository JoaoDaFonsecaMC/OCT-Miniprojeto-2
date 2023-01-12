function [c, ceq] = confun1(x)
% Nonlinear inequality constraints:
c = [-8+x(1)-x(2)+x(3)-x(4)+x(1)^(2)+x(2)^(2)+x(3)^(2)+x(4)^(2); 
 -10-x(1)-x(4)+x(1)^(2)+2*x(2)^(2)+x(3)^(2)+2*x(4)^(2)];
%  nonlinear equality constraints:
ceq = [5-2*x(1)+x(2)+x(4)-2*x(1)^(2)-x(2)^(2)-x(3)^(2)];
end