%main script task 3 

%Matrixes construction
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
mu_q=1; % between 1 and 100 including 1 e 100
mu_r=0.1; % between 0 and 1 excluding 0 e 1
Q=mu_q*eye(6);
R=mu_r*eye(3);

%LQR algorithm
[K,P,E]=lqr(A,B,Q,R); % you need Control System Toolbox to run it 

%Waypoints coordinates and times

Way7= [39 49 25.71 ; 7 29 35];
Way8= [39 50 34.82 ; 7 29 37];
Way9= [39 51 33.38 ; 7 29 39];
Way10=[39 52 39.86 ; 7 29 41];
H7=400;
H8=500;
H9=600;
H10=600;
t7=0*3600; %times from hours to seconds
t8=0.035*3600;
t9=0.070*3600;
t10=0.080*3600;

%Waypoints in degress
Way7_deg=dms2degrees(Way7); % you need Mapping Toolbox to run it 
Way8_deg=dms2degrees(Way8);
Way9_deg=dms2degrees(Way9);
Way10_deg=dms2degrees(Way10);

%Transform West with minus symbol
Way7_deg(2)=-Way7_deg(2);
Way8_deg(2)=-Way8_deg(2);
Way9_deg(2)=-Way9_deg(2);
Way10_deg(2)=-Way10_deg(2);

%Waypoints in geocentric coordinates 
[X7,Y7,Z7]= geodetic_to_geocentric (Way7_deg(2),Way7_deg(1),H7);
[X8,Y8,Z8]= geodetic_to_geocentric (Way8_deg(2),Way8_deg(1),H8);
[X9,Y9,Z9]= geodetic_to_geocentric (Way9_deg(2),Way9_deg(1),H9);
[X10,Y10,Z10]= geodetic_to_geocentric (Way10_deg(2),Way10_deg(1),H10);

%Linear regression to build Yref
t=t7;
step=0.01;
y_reftotal=[];
Matrix1=[A B; C zeros(3)];
Matrix2=[zeros(6,3); eye(3)];
Asterix_matrix=[];
j=1;

while(t<t8)
 Y_refnew7= [X7+(X8-X7)/(t8-t7)*(t-t7);
             Y7+(Y8-Y7)/(t8-t7)*(t-t7);
             Z7+(Z8-Z7)/(t8-t7)*(t-t7)];
 y_reftotal=cat(2,y_reftotal,Y_refnew7); % cat concanates arrays
 New_asterix=inv(Matrix1)*Matrix2*y_reftotal(:,j);
 Asterix_matrix=cat(2,Asterix_matrix,New_asterix);
 t=t+step;
 j=j+1;
end

while(t<t9)
 Y_refnew8=[X8+(X9-X8)/(t9-t8)*(t-t8);
             Y8+(Y9-Y8)/(t9-t8)*(t-t8);
             Z8+(Z9-Z8)/(t9-t8)*(t-t8)];
 y_reftotal=cat(2,y_reftotal,Y_refnew8); % cat concanates arrays
 New_asterix=inv(Matrix1)*Matrix2*y_reftotal(:,j);
 Asterix_matrix=cat(2,Asterix_matrix,New_asterix);
 t=t+step;
 j=j+1;
end

while(t<t10)
    Y_refnew9=[X9+(X10-X9)/(t10-t9)*(t-t9);
             Y9+(Y10-Y9)/(t10-t9)*(t-t9);
             Z9+(Z10-Z9)/(t10-t9)*(t-t9)];
 y_reftotal=cat(2,y_reftotal,Y_refnew9); % cat concanates arrays
 New_asterix=inv(Matrix1)*Matrix2*y_reftotal(:,j);
 Asterix_matrix=cat(2,Asterix_matrix,New_asterix);
 t=t+step;
 j=j+1;
end

mu_w=0.2; % value between 0 and 1
w=mu_w*eye(3);

Ad=expm(A*step); % matrix exponential
N=2;
Bd=(eye(6)*step-A^(N-1)*step^N)/factorial(N)*B;
v_0=[0 0 0 ]'; 
x_ast=Asterix_matrix(1:6,:);
uc_ast=Asterix_matrix(7:9,:);
k=1;
xtotal=[];
uctotal=[];
x0=cat(1,y_reftotal(:,1),v_0);
xk=x0;

t=t7;
while(t<=t8)
 uc_k=inv(w)*(w*uc_ast(:,k)-B'*P*(xk-x_ast(:,k)));
 xk=Ad*xk+Bd*uc_k;
 k=k+1;
 t=t+step;
 xtotal=cat(2,xtotal,xk);
 uctotal=cat(2,uctotal,uc_k);
end

while(t<=t9)
 uc_k=inv(w)*(w*uc_ast(:,k)-B'*P*(xk-x_ast(:,k)));
 xk=Ad*xk+Bd*uc_k;
 k=k+1;
 t=t+step;
  xtotal=cat(2,xtotal,xk);
 uctotal=cat(2,uctotal,uc_k);
end

while(t<=t10)
 uc_k=inv(w)*(w*uc_ast(:,k)-B'*P*(xk-x_ast(:,k)));
 xk=Ad*xk+Bd*uc_k;
 k=k+1;
 t=t+step;
  xtotal=cat(2,xtotal,xk);
 uctotal=cat(2,uctotal,uc_k);
end

time=t7:step:t10;
%Plots
figure ()
plot3(xtotal(1,:) , xtotal(2,:) , xtotal(3,:));
xlabel('X[Km]')
ylabel('Y[Km]')
zlabel('Z[Km]')
grid on
hold on
plot3(X7 , Y7 , Z7, 'o');
plot3(X8 , Y8 , Z8, 'o');
plot3(X9 , Y9 , Z9, 'o');
plot3(X10 , Y10 , Z10, 'o');

figure ()
plot(time,xtotal(4,:))
title('vx in time')

figure ()
plot(time,xtotal(5,:))
title('vy in time')

figure ()
plot(time,xtotal(6,:))
title('vz in time')

figure ()
plot(time,uctotal(1,:))
title('ax in time')

figure ()
plot(time,uctotal(2,:))
title('ay in time')

figure ()
plot(time,uctotal(3,:))
title('az in time')


