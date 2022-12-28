function [u_dot] = f_u(phi,v,r,miu,T,m)

    u_dot = v^(2)/r - miu/r^(2) + T * sin(phi)/m
end