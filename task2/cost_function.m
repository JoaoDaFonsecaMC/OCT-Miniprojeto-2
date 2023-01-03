function [J] = cost_function (phi)
%parameters
N = 100;
t_0 = 0;
t_f = 1;

%time step
h = (t_f-t_0)/N;

%initial values
u_0 = 0;
v_0 = 1/sqrt(1.1);
r_0 = 1.1;
theta_0 = 0.1;

U = zeros(N+1);
V = zeros(N+1);
R = zeros(N+1);
Theta = zeros(N+1);

%initialization
U(1) = u_0;
V(1) = v_0;
R(1) = r_0;
Theta(1) = theta_0;

%cosntans
m0 = 1;
m_dot = -0.07487;
T = 0.1405;
miu = 1;

consts = [m0, m_dot, T, miu];

%discretization
[U,V,R,Theta]=euler_implicit_fix_point(phi,@f_u,@f_v,@f_r,@f_theta,@mass,u_0,v_0,r_0,theta_0,t_f,h,N,consts);

%cost function to optimize
J = [(-1)*(0.5*(U(N+1)^(2)+V(N+1)^(2))-1/R(N+1))];

end