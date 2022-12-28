function [u_f,v_f,r_f]=euler_implicit_fix_point(phi,u_dot,v_dot,r_dot,f_m,u_0,v_0,r_0,t_max,consts)

delta_t = 0.01;

m0 = consts(1);
m_dot = consts(2);
T = consts(3);
miu = consts(4);

% [t_h,u_h,vett_it_newton]=eulero_indietro_pto_fisso (f, t_max, valore_iniziale, delta_t)
%
% Metodo di Eulero all'indietro per la soluzione di equazioni differenziali
% ordinarie con metodo di punto fisso per il calcolo della soluzione ad ogni
% step
%
% Parametri di ingresso:
%
% f                 funzione di t, u associata al problema di Cauchy
% t_max             estremo di integrazione
% valore_iniziale   valore della soluzione al tempo iniziale t0 = 0
% delta_t           passo di integrazione
% 
%
% Parametri di uscita:
%
% t_h               vettore dei tempi in cui la soluzione viene calcolata
% u_h               vettore delle soluzioni calcolate in t_h
% vett_it_pf        vettore delle iterazioni del metodo di punto fisso 
%                   ad ogni passo

t0=0;

t_h=t0:delta_t:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=length(t_h);

u_h=zeros(1,N_istanti);
v_h=zeros(1,N_istanti);
r_h=zeros(1,N_istanti);

% ciclo iterativo che calcola u_(n+1)=u_n+h*f_(n+1) . Ad ogni iterazione temporale devo eseguire delle sottoiterazioni
% di punto fisso per il calcolo di u_(n+1): u_(n+1)^(k+1) = u_n + h * f_( t_(n+1) , u_(n+1)^k ). Per garantire la
% convergenza di tale metodo la derivata della funzione di iterazione phi(x) u_n + h * f_( t_(n+1) , x ) deve avere
% derivata minore di 1 in modulo. Questo introdurra' delle condizioni su h.

u_h(1)=u_0;
v_h(1)=v_0;
r_h(1)=r_0;

% parametri per le iterazioni di punto fisso
N_max=100;
toll=1e-5;

vett_it_pf=zeros(1,N_istanti);

for it=2:N_istanti
    
    % preparo le variabili per le sottoiterazioni
    
    u_old=u_h(it-1);
    v_old=v_h(it-1);
    r_old=r_h(it-1);
    
    t_pf=t_h(it);
    
    m = f_m(t_pf,m0,m_dot);
        
    psi_u=@(u,phi,v,r,miu,T,m) u_old + delta_t * u_dot( phi,v,r,miu,T,m );
    psi_v=@(v,phi,u,r,T,m) v_old + delta_t * v_dot( phi,v,u,r,T,m );
    psi_r=@(r,u) r_old + delta_t * r_dot( u );
    
    % sottoiterazioni
    
    [u_pf, it_pf] = ptofis_u(u_old, psi_u, N_max, toll,phi,v_old,r_old,miu,T,m);
    [v_pf, it_pf] = ptofis_v(v_old, psi_v, N_max, toll,phi,u_old,r_old,T,m);
    [r_pf, it_pf] = ptofis_r(r_old, psi_r, N_max, toll,u_old);
    
    u_h(it)=u_pf(end);
    v_h(it)=v_pf(end);
    r_h(it)=r_pf(end);
    
    % tengo traccia dei valori di it_pf per valutare la convergenza delle iterazioni di punto fisso
    vett_it_pf(it)=it_pf;
    
end
u_f = u_h(end);
v_f = v_h(end);
r_f = r_h(end);
end