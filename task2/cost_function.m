function [J] = cost_function (tf,u,v,r)

J = (-1)*(0.5*(u(tf)^(2)+v(tf)^(2))-1/(r(tf)))

end