function [Euler_diff] = euler_implicit(Phi,U,V,R,Theta,consts)
    %load constants
    N = consts(1);
    h = consts(2);
    T = consts(3);
    miu = consts(4);
    m0 = consts(5);
    m_dot = consts(6);
    
    k1 = zeros(4);
    k2 = zeros(4);
    
    Euler_diff = zeros(N,4);
    
    %main cycle
     for i=1:N
         %time
         t = h*i;
         
         %u
         [k1(1)] = f_u(Phi(i),V(i),R(i),miu,T,mass(t,m0,m_dot));
         [k2(1)] = f_u(Phi(i+1),V(i+1),R(i+1),miu,T,mass(t+h,m0,m_dot));
        
         Euler_diff(i,1) = U(i+1) -( U(i) + (h/2)*(k1(1)+k2(1)));
         
         %v
         [k1(2)] = f_v(Phi(i),V(i),U(i),R(i),T,mass(t,m0,m_dot));
         [k2(2)] = f_v(Phi(i+1),V(i+1),U(i+1),R(i+1),T,mass(t+h,m0,m_dot));
        
         Euler_diff(i,2) = V(i+1) -( V(i) + (h/2)*(k1(2)+k2(2)));
         
         %r
         [k1(3)] =  f_r(U(i));
         [k2(3)] = f_r(U(i+1));
        
         Euler_diff(i,3) = R(i+1) -( R(i) + (h/2)*(k1(3)+k2(3)));
         
         %theta
         [k1(4)] = f_theta(V(i),R(i));
         [k2(4)] = f_theta(V(i+1),R(i+1));
        
         Euler_diff(i,4) = Theta(i+1) -( Theta(i) + (h/2)*(k1(4)+k2(4)));
     end
end