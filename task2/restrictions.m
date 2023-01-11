function [c, ceq ] = restrictions(x)
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
Phi = x(:,1);
U = x(:,2);
V = x(:,3);
R = x(:,4);
Theta = x(:,5);

X1 = [U(1) V(1) R(1) Theta(1)];
x0 = [u_0 v_0 r_0 theta_0];

[Euler_diff] = euler_implicit(Phi,U,V,R,Theta,consts);

% Nonlinear inequality constraints:
c = [Phi - 0.5; -Phi - 0.5];
% no nonlinear equality constraints:
ceq = [X1 - x0; Euler_diff];

end