function [J] = cost_function (x)
%loading values
Phi = x(:,1);
U = x(:,2);
V = x(:,3);
R = x(:,4);
Theta = x(:,5);

%cost function to optimize
J = [(-1)*(0.5*(U(end)^(2)+V(end)^(2))-1/R(end))];

end