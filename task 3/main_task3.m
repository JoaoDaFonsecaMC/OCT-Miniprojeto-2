format long
A= [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];
B= [0 0 0;
    0 0 0;
    0 0 0;
    1 0 0;
    0 1 0;
    0 0 1];
C= [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0];
MewQ=1; % between 1 and 100 including 1 e 100
MewR=0.1; % between 0 and 1 excluding 0 e 1
Q=MewQ*eye(6);
R=MewR*eye(3);
[K,P,E]=lqr(A,B,Q,R) % you need Control System Toolbox to run it 
Way7= [39 59 15.04 ; 7 29 45];
Way8= [40 01 17.12 ; 7 29 47];
Way9= [40 03 45.92 ; 7 29 49];
Way10=[40 05 31.38 ; 7 29 51];
H7=800;
H8=800;
H9=800;
H10=800;
t7=0.210*3600;
t8=0.245*3600;
t9=0.280*3600;
t10=0.325*3600;
