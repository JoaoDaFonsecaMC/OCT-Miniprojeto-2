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

Way1= [39 49 25.71 ; 7 29 35];
Way2= [39 50 34.82 ; 7 29 37];
Way3= [39 51 33.38 ; 7 29 39];
Way4=[39 52 39.86 ; 7 29 41];
H1=400;
H2=500;
H3=600;
H4=600;
t1=0*3600; %times from hours to seconds
t2=0.035*3600;
t3=0.070*3600;
t4=0.080*3600;

%Waypoints in degress
Way1_deg=dms2degrees(Way1); % you need Mapping Toolbox to run it 
Way2_deg=dms2degrees(Way2);
Way3_deg=dms2degrees(Way3);
Way4_deg=dms2degrees(Way4);

%Transform West with minus symbol
Way1_deg(2)=-Way1_deg(2);
Way2_deg(2)=-Way2_deg(2);
Way3_deg(2)=-Way3_deg(2);
Way4_deg(2)=-Way4_deg(2);

%Waypoints in geocentric coordinates 
[X1,Y1,Z1]= geodetic_to_geocentric (Way1_deg(2),Way1_deg(1),H1);
[X2,Y2,Z2]= geodetic_to_geocentric (Way2_deg(2),Way2_deg(1),H2);
[X3,Y3,Z3]= geodetic_to_geocentric (Way3_deg(2),Way3_deg(1),H3);
[X4,Y4,Z4]= geodetic_to_geocentric (Way4_deg(2),Way4_deg(1),H4);

%Linear regression to build Yref
t=t1;
step=0.01;
y_reftotal=[];
Matrix1=[A B; C zeros(3)];
Matrix2=[zeros(6,3); eye(3)];
Asterix_matrix=[];
j=1;

while(t<t2)
 Y_refnew1= [X1+(X2-X1)/(t2-t1)*(t-t1);
             Y1+(Y2-Y1)/(t2-t1)*(t-t1);
             Z1+(Z2-Z1)/(t2-t1)*(t-t1)];
 y_reftotal=cat(2,y_reftotal,Y_refnew1); % cat concanates arrays
 New_asterix=inv(Matrix1)*Matrix2*y_reftotal(:,j);
 Asterix_matrix=cat(2,Asterix_matrix,New_asterix);
 t=t+step;
 j=j+1;
end

while(t<t3)
 Y_refnew2=[X2+(X3-X2)/(t3-t2)*(t-t2);
             Y2+(Y3-Y2)/(t3-t2)*(t-t2);
             Z2+(Z3-Z2)/(t3-t2)*(t-t2)];
 y_reftotal=cat(2,y_reftotal,Y_refnew2); % cat concanates arrays
 New_asterix=inv(Matrix1)*Matrix2*y_reftotal(:,j);
 Asterix_matrix=cat(2,Asterix_matrix,New_asterix);
 t=t+step;
 j=j+1;
end

while(t<t4)
    Y_refnew3=[X3+(X4-X3)/(t4-t3)*(t-t3);
             Y3+(Y4-Y3)/(t4-t3)*(t-t3);
             Z3+(Z4-Z3)/(t4-t3)*(t-t3)];
 y_reftotal=cat(2,y_reftotal,Y_refnew3); % cat concanates arrays
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

t=t1;
while(t<=t2)
 uc_k=inv(w)*(w*uc_ast(:,k)-B'*P*(xk-x_ast(:,k)));
 xk=Ad*xk+Bd*uc_k;
 k=k+1;
 t=t+step;
 xtotal=cat(2,xtotal,xk);
 uctotal=cat(2,uctotal,uc_k);
end

while(t<=t3)
 uc_k=inv(w)*(w*uc_ast(:,k)-B'*P*(xk-x_ast(:,k)));
 xk=Ad*xk+Bd*uc_k;
 k=k+1;
 t=t+step;
  xtotal=cat(2,xtotal,xk);
 uctotal=cat(2,uctotal,uc_k);
end

while(t<=t4)
 uc_k=inv(w)*(w*uc_ast(:,k)-B'*P*(xk-x_ast(:,k)));
 xk=Ad*xk+Bd*uc_k;
 k=k+1;
 t=t+step;
  xtotal=cat(2,xtotal,xk);
 uctotal=cat(2,uctotal,uc_k);
end

time=t1:step:t4;
%Plots
figure ()
plot3(xtotal(1,:) , xtotal(2,:) , xtotal(3,:));
xlabel('X[Km]')
ylabel('Y[Km]')
zlabel('Z[Km]')
title('3D Navigation')
grid on
hold on
plot3(X1 , Y1 , Z1, 'o');
plot3(X2 , Y2 , Z2, 'o');
plot3(X3 , Y3 , Z3, 'o');
plot3(X4 , Y4 , Z4, 'o');

figure ()
plot(time,xtotal(4,:))
xlabel ('time [s]')
ylabel ('vx [km/s]')
title('vx in time')

figure ()
plot(time,xtotal(5,:))
xlabel ('time [s]')
ylabel ('vy [km/s]')
title('vy in time')

figure ()
plot(time,xtotal(6,:))
xlabel ('time [s]')
ylabel ('vz [km/s]')
title('vz in time')

figure ()
plot(time,uctotal(1,:))
xlabel ('time [s]')
ylabel ('ax [km/s^2]')
title('ax in time')

figure ()
plot(time,uctotal(2,:))
xlabel ('time [s]')
ylabel ('ay [km/s^2]')
title('ay in time')

figure ()
plot(time,uctotal(3,:))
xlabel ('time [s]')
ylabel ('az [km/s^2]')
title('az in time')


