%% Magnetically steerable catheter magnetic localization simulation
% model: ideal dipole model
% sensor configuration: 1. rectangular sensor array
%                       2. circular sensor array
% variables: distance 

% Start with 4x4 arrays in rectangular sensor plate array
clear;close all;clc;

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

OD_mag = 3.175e-3;  % outer diameter of PM  [m]
ID_mag = OD_mag/2;  % inner diameter of PM  [m]
length_mag = OD_mag;   % length of PM  [m]
V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
Br = 1.2; % magnetic remanance

M = V_pm*Br/mu_0;

% Whole procedure: set the PM location and orientation in advance, compute the
% analytical magnetic field in all 16 sensors locations in hall effect sensor measuring direction. 
% Use the 16 sensor measurements + manual added noise with 
% different level of magnitude and perform L-M least square optimization
% algorithm to estimate the location and orientation

% Assumption: stationary PM magnet (later with trajectory in space over time and estimate the trajectory)
%             ideal dipole model

% Sensor array parameters
% All coordinates here
SID = 0.01;  % sensor interval distance [m]
% SHD = 0.01;  % sensor array height w.r.t working range of PM [m]

% 4x4 Square sensor array confuguration
n_sens = 16;
x_s = [-1.5*SID   -0.5*SID   0.5*SID   1.5*SID   -1.5*SID   -0.5*SID   0.5*SID   1.5*SID  ...
    -1.5*SID   -0.5*SID   0.5*SID   1.5*SID   -1.5*SID   -0.5*SID   0.5*SID   1.5*SID ];
y_s = [1.5*SID   1.5*SID   1.5*SID   1.5*SID   0.5*SID   0.5*SID   0.5*SID   0.5*SID  ...
    -0.5*SID   -0.5*SID   -0.5*SID   -0.5*SID   -1.5*SID   -1.5*SID   -1.5*SID   -1.5*SID ];
z_s = zeros(1,n_sens);
xyz_s = [x_s;y_s;z_s];  % i-th column is i-th sensor location (x_i,y_i,z_i)

% define measurement direction vector (3 x m)
% i-th column shows the [x,y,z] components of the measurement at i-th
% sensor
meas_dir_s = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ;
              0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ;
              1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 ];
meas_dir_s = meas_dir_s./sqrt(sum(meas_dir_s.^2));  % normalization to unit norm in each direction


% Preset PM dipole location and orientation variable
% (x0,y0,z0)  (sin(theta0)cos(phi0) , sin(theta0)sin(phi0) , cos(theta0))
x0 = 0;  % [m]
y0 = 0;  % [m]
z0 = 0.1;  % [m]
theta0 = 0;  % [rad]
phi0 = 0;  % [rad]
% theta0 = pi/2;  % [rad]
% phi0 = pi/2;  % [rad]
xyz_pm = [x0;y0;z0];
pm_all = [x0;y0;z0;theta0;phi0];

figure
plot3(x_s,y_s,z_s,'rx')
hold on 
% plot3([x0,x0+1e-2*sin(theta0)*cos(phi0)],[y0,y0+1e-2*sin(theta0)*sin(phi0)],[z0,z0+1e-2*cos(theta0)],'k')
% quiver3(x0,y0,z0,[1e-2*sin(theta0)*cos(phi0)],[1e-2*sin(theta0)*sin(phi0)],[1e-2*cos(theta0)],'LineWidth',2)

% Disturbance field G
G_x = 0;
G_y = 0;
G_z = 0;
G = [G_x;G_y;G_z];

% Given known PM position and orientation, magnetic moment, sensor array
% configuration, compute the magnetic field at all sensor locations
% B_sx/B_sy/B_sz are 1x16 vector with i-th element being i-th sensor x/y/z field
B_measure_scalar = PM_forward_field(M,G,pm_all,xyz_s,meas_dir_s);

% Timeline 0~t_final
% first assume stationary PM s.t. at each moment the forward generated
% field is constant
t_final = 10;  % [s]
dt = 0.1;  % sampling time [s]
t = 0:dt:t_final;  % time history
n_t = length(t);  % length of time data

% Simulate measurement over time
% Here consider square array s.t. measurement is B_sz
% B_meas_ideal rows: time index, columns: i-th sensor
B_meas_ideal = ones(n_t,1)*B_measure_scalar;

% Add noise into measurement on purpose with different level of magnitude
noise_mu = 0;  % mean
noise_sigma = 1e-5;  % variance
noise = noise_mu + noise_sigma*randn(n_t,n_sens);  % generating random noise according to given distribution

% Final corrupted measurement 
B_meas = B_meas_ideal + noise;

% some plots
figure
plot(t,B_meas(:,6))
hold on 
plot([0,t(end)],[B_meas_ideal(1,6) B_meas_ideal(1,6)])
xlabel('Time[s]')
ylabel('Sensor measurement')
legend('corrupted','ideal')

figure
plot(t,B_meas_ideal)

%% Backward estimation phase
% Given measurement of all sensors over a time horizon
% consider at one time, the measurement B_meas is 1 x n_sens vector

% Objective function: f_err = || B_est - B_meas ||^2
% Variable: pm_hat = (x_hat,y_hat,z_hat,theta_hat,phi_hat) as estimation of the position and orientation

% initial guess on parameter estimation at time 0
pm_hat_0 = [0.01;0.01;0.00;0.2;0.2];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0)

% Compute Analytic Jacobian matrix： J <- R^(m x n) here
J = Analytic_Jacobian_Mag_sqr_array(M,pm_hat_0,xyz_s);


% Levenberg–Marquardt least square algorithm solving for optimal position/orientation estimation
% test
B_meas_test = B_measure_scalar';
pm_hat_opt = PM_backward_estimation(B_meas_test,pm_hat_0,M,G,xyz_s,meas_dir_s); 

% weighting among parameters to be added!

%% Timeline estimation
% use the current estimation for next step initial guess
% record the requiered computation time
pm_est_history = zeros(n_t,5);
pm_err_history = zeros(n_t,5);
comp_time_history = zeros(n_t,1);
pm_last = pm_hat_0;
for i = 1:n_t
    tic
    B_meas_i = B_meas(i,:)';
    pm_last = PM_backward_estimation(B_meas_i,pm_last,M,G,xyz_s,meas_dir_s); 
    pm_est_history(i,:) = pm_last;
    pm_err_history(i,:) = pm_last - pm_all; 
    comp_time_history(i) = toc;
end


% plot the estimation over time
figure()
plot(t,pm_err_history(:,1))
xlabel('Time [s]')
ylabel(' X position [m]')
title('X position estimation error over time')

figure()
plot(t,pm_err_history(:,2))
xlabel('Time [s]')
ylabel(' Y position [m]')
title('Y position estimation error over time')

figure()
plot(t,pm_err_history(:,3))
xlabel('Time [s]')
ylabel(' Z position [m]')
title('Z position estimation error over time')

figure()
plot(t,pm_err_history(:,4))
xlabel('Time [s]')
ylabel(' Theta orientation [rad]')
title('Theta orientation estimation error over time')

figure()
plot(t,pm_err_history(:,5))
xlabel('Time [s]')
ylabel(' Phi orientation [rad]')
title('Phi orientation estimation error over time')


figure
plot(t,comp_time_history)
xlabel('Time [s]')
ylabel('Computation time [s]')

%% Error analysis
% Test Sensor Height Distance v.s. position/orientation error




% Test signal/noise ratio on purpose




%% random test 
y_meas = B_sz';
pm_init = pm_hat_0;





