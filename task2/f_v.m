function [v_dot] = f_v(phi,v,u,r,T,m)

    v_dot = -u*v/r + T * cos(phi)/m
end