function [J] = cost_function (phi)

%max time
t_max = 1;

%initial values
u_0 = 0;
v_0 = 1/sqrt(1.1);
r_0 = 1.1;

%cosntans
m0 = 1;
m_dot = -0.07487;
T = 0.1405;
miu = 1;

consts = [m0, m_dot, T, miu];

%discretization
[u_h,v_h,r_h]=euler_implicit_fix_point(phi,@f_u,@f_v,@f_r,@mass,u_0,v_0,r_0,t_max,consts);

%cost function to optimize
J = [(-1)*(0.5*(u_h^(2)+v_h^(2))-1/r_h)];

end