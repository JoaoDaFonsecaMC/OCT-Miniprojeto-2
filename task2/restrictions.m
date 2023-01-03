function [c, ceq ] = restrictions(phi)

% Nonlinear inequality constraints:
c = [phi(1,:) - 0.5; -phi(1,:) - 0.5];
% no nonlinear equality constraints:
ceq = [];

end