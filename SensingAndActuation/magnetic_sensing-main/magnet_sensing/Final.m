clear all
clc

%setting variables
g = 9.81; %m/s^2
L = 11.2e-2; %m
mp = 381e-3; %kg
Ip = 0.00616; %kg m2
mw = 36e-3; %kg
i_cmw = 7.46e-6; %kg m^2
rw = 2.1e-2; %m
R = 4.4; %ohm
kb = 0.495; %Vs/rad
kt = 0.470; %Nm/A

q1 = L*mp;
q2 = Ip + (L^2)*mp;
q3 = (kb*kt)/R;
q4 = -mp - mw - i_cmw/rw^2;

%organizing matrices
A12 = (g*q1*q4)/(q1^2+q2*q4);
A22 = (q3 * (q1- q4*rw))/((q1^2+q2*q4)*rw);
A42 = (q3*(-q1 + q4*rw))/((q1^2+q2*q4)*(rw^2));
A14 = (g*q1^2)/(q1^2+q2*q4);
A24 = -(q3 * (q2+ q1*rw))/((q1^2+q2*q4)*rw);
A44 = (q3*(q2 + q1*rw))/((q1^2+q2*q4)*(rw^2));
B2 = (kt*(q1-q4*rw))/(R*(q1^2+q2*q4)*rw);
B4 = 20.6;
% B4 = (kt*(-q2+q1*rw))/(R*(q1^2-q2*q4)*rw);

A = [0, 1, 0, 0; 
    A12, A22, 0 , A42;
    0, 0, 0, 1;
    A14, A24, 0 , A44]

B = [0;B2;0;B4]
C = [1,1,1,1]
D = 0;

I = [1, 0, 0, 0;
     0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];
sysu = ss(A,B,C,D);
syms s lambda
CE = det(lambda*I - A) %characteristic equation 
%Controllability and observability
Cont = ctrb(sysu) 
controllability = rank(Cont);
if controllability == 4
    disp('System is controllable')
else
    disp('System is NOT controllable')
end
Ob = obsv(sysu); 
observability = rank(Ob);
if observability == 4
    disp('System is observable')
el
end


%Root Locus plotting
[tf_num, tf_den] = ss2tf(A,B,C,D)
sys = tf(tf_num, tf_den)
figure(1), rlocus(sys)
sys_cl = feedback(sys,1)
figure(1), step(sys_cl)
grid on

% % specify desired zeta and wn (iterate as needed)
% zeta = 0.8; wn = 0.25;
% 
% Ds = expand((s*s + 2*zeta*wn*s + wn*wn)*(s+zeta*wn));
% ps = coeffs(Ds,s);
% pn = double(ps);
% p = roots(pn);
% K = place(A,B,p);
% sysc = ss((A-B*K),B,C,D);
% 
% 
% 
% % uncontrolled system
% sysu = ss(A,B,C,D);
% %state variables x = [alpha, alpha_dot, x, x_dot]
% initial_condition = [deg2rad(20), 0, 10, 0];
% %input vector = 
% % input_vector = [0; 12; 0; 12];
% [yu,tu,xu] = initial(sysu,initial_condition,5);
% plot(tu,xu)
