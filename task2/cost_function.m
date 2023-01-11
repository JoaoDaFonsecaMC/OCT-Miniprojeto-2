function [J] = cost_function (X)
%loading values
Phi = X(1,:);
U = X(2,:)
V = X(3,:)
R = X(4,:)
Theta = X(5,:)

%cost function to optimize
J = [(-1)*(0.5*(U(end)^(2)+V(end)^(2))-1/R(end))];

end