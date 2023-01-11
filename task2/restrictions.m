function [c, ceq ] = restrictions(X)
%parameters
N = 100;
t_0 = 0;
t_f = 1;

%initial values
u_0 = 0;
v_0 = 1/sqrt(1.1);
r_0 = 1.1;
theta_0 = 0.1;
m0 = 1;
m_dot = -0.07487;
T = 0.1405;
miu = 1;

%time step
h = (t_f-t_0)/N;

consts = [N h T miu m0 m_dot];

%loading values
Phi = X(1,:);
U = X(2,:)
V = X(3,:)
R = X(4,:)
Theta = X(5,:)

[Euler_diff] = euler_implicit(Phi,U,V,R,Theta,consts)

% Nonlinear inequality constraints:
c = [Phi(1,:) - 0.5; -Phi(1,:) - 0.5];
% no nonlinear equality constraints:
ceq = [X(2,1) - u_0; X(3,1) - v_0; X(4,1) - r_0; X(5,1) - theta_0; Euler_diff];

end