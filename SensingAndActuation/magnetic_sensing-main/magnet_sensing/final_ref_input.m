% FSFB with reference input for three-state (hydraulic) system
% rgl 11-27-22
% no disturbance input
clear all
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
    A14, A24, 0 , A44];

B = [0;B2;0;B4];

C = [1, 0, 1, 0];
% C = [1, 0, 0, 0;
%     0, 1, 0, 0;
%     0, 0, 1, 0;
%     0, 0, 0, 1];

D = 0;



% % for adding in disturbance into tank 1
% Bd = [1;0;0];
% uncontrolled system
sysu = ss(A,B,C,D);
% set poles
syms s 
% specify desired zeta and wn (iterate as needed)
Q = [1,0,0,0;
     0,0.001,0,0;
     0,0,100,0;
     0,0,0,1]; 
R = 1;
[K_LQR,S_LQR,P_LQR] = lqr(A,B,Q,R);
K = place(A,B,P_LQR);
% compare response of OL system and state-feedback regulator
Ac = A-B*K;
Bc = [0;0;0;0];

sysc = ss(Ac,B,C,D);

Nbar = rscale(A,B,C,D,K);
% set initial reference value
ro = deg2rad(20);
% solve for initial conditions for given reference value
% (from equilibrium of the closed-loop system equations)
xo = [ro, 0, ro, 0]';
% xo = [deg2rad(20),0,deg2rad(20),0]';
% note, 1 gallon = 0.003785 m^3
% time array
t = 0:0.001:1;
% compute the desired reference values over time
% note: these need to be realistic and could cause 
% system states or control input to
% take values that are not realistic
T = 0.4;
tstart = 0.1;
for i = 1:length(t)
    if t(i)>tstart
        r(i) = xo(1)*(1+square(2*pi*(t(i)-tstart)/T));
    else
        r(i) = xo(1);
    end
end

% compute the reference for lsim() simulation
% NOTE: This ur will multipled by B defined above
ur = Nbar*r;
uro = Nbar*ro;
uco = -K*xo;
% simulate the controlled system with reference input
% and specified initial conditions xo
[yc,tc,xc] = lsim(sysc,ur,t,xo);
% compute the control inputs
uc = -K*xc';
ucm = max(uc);
sprintf('ucm = %4.2f',ucm)
u = uc + ur;
%
figure(1)
subplot(2,1,1), plot(tc,r,tc,xc(:,1),tc,xc(:,3))
legend('r', 'alpha', 'x')
%ylim([0.04,0.06])
% subplot(3,1,2), plot(tc,yc(:,2))
% legend('Voltage')
subplot(2,1,2), plot(t,uc,t,-ur), legend('u_c','-u_r')
%ylim([5,6.5])