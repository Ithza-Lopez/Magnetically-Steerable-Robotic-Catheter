%% Magnetically steerable catheter magnetic localization simulation
% model: ideal dipole model
% sensor configuration: 1. rectangular sensor array
%                       2. circular sensor array
% variables: distance 

% Start with 4x4 arrays in circular sensor plate array
clear;close all;clc;

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

%smallmagnet
OD_mag = 3.175e-3;  % outer diameter of PM  [m]
ID_mag = OD_mag/2;  % inner diameter of PM  [m]
length_mag = OD_mag*2;   % length of PM  [m]
V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM

% %cube
% V_pm = (0.25*0.0254)^3;

Br = 1.48; % magnetic remanance

M = V_pm*Br/mu_0;  % estimated magnetic moment


% Whole procedure: set the PM location and orientation in advance, compute the
% analytical magnetic field in all 16 sensors locations in hall effect sensor measuring direction. 
% Use the 16 sensor measurements + manual added noise with 
% different level of magnitude and perform L-M least square optimization
% algorithm to estimate the location and orientation

% Assumption: stationary PM magnet (later with trajectory in space over time and estimate the trajectory)
%             ideal dipole model

% Sensor array parameters
% All coordinates here
L = 0.025;  % sensor interval distance [m]
R = (79.67/2)*1e-3;  % circular sensor array radius [m]

%% Loading measured data

filename = '1-28-23 Small Magnet Centered 2G 2';
B_meas = B_array(filename);

%% 4x4 Circular sensor array confuguration
n_sens = 16;

magnet_distance = 41.39; %mm
dim_repeat = ones(1,3)*magnet_distance; 

%all sensors working
x_s = [repelem((zeros(1,3)),4), repelem((-dim_repeat),4), repelem((dim_repeat),4), repelem((zeros(1,3)),4)]*1e-3; % 6 zeros and alternating x pos and neg
y_s = [repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3),repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3),repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3),repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3)]*1e-3; % decreasing from 126.4 to 51.4 and increaseing to 126.4.
z_s = [repelem((dim_repeat),4),repelem((zeros(1,3)),8), repelem((-dim_repeat),4)]*1e-3; %alternating z pos and neg, and zeros
xyz_s = [x_s;y_s;z_s];  % i-th column is i-th sensor location (x_i,y_i,z_i)

% define measurement direction vector (3 x m)
% i-th column shows the [x,y,z] components of the measurement at i-th
% sensor
%all sensor working
meas_dir_s = [eye(3).*[-1;-1;1],...
              eye(3).*[-1;-1;1],...
              eye(3).*[-1;-1;1],...
              eye(3).*[-1;-1;1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[1;-1;1],...
              fliplr(eye(3)).*[1;-1;1],...
              fliplr(eye(3)).*[1;-1;1],...
              fliplr(eye(3)).*[1;-1;1],...
              eye(3).*[1;-1;-1],...
              eye(3).*[1;-1;-1],...
              eye(3).*[1;-1;-1],...
              eye(3).*[1;-1;-1]];

meas_dir_s = meas_dir_s./sqrt(sum(meas_dir_s.^2));  % normalization to unit norm in each direction


%% Preset PM dipole location and orientation variable
% (x0,y0,z0)  (sin(theta0)cos(phi0) , sin(theta0)sin(phi0) , cos(theta0))
x0 = 0.0;  % [m]
y0 = 0.089;  % [m]
z0 = 0.0;  % [m]
theta0 = pi/2;  % [rad]
phi0 = pi/2;  % [rad]
% theta0 = pi/2;  % [rad]
% phi0 = pi/2;  % [rad]
xyz_pm = [x0;y0;z0];
m0 = M;
pm_all = [x0;y0;z0;theta0;phi0;m0];
% pm_all = [0.00163458779397101;0.0860972360316718;0.000904777777496003;1.70719589598097;1.55816453265289];

% Disturbance field G
G_x = 0;
G_y = 0;
G_z = 0;
G = [G_x;G_y;G_z];

% saturation obtained from actual Hall Effect sensor
sat = 0.1;  % [T]

% Given known PM position and orientation, magnetic moment, sensor array
% configuration, compute the magnetic field at all sensor locations
% B_sx/B_sy/B_sz are 1x16 vector with i-th element being i-th sensor x/y/z field
% compute the actual saclar composite measurement (1 x m) at all sensor
% locations
B_measure_theoretical = PM_forward_field(G,pm_all,xyz_s,meas_dir_s);

% %% Timeline 0~t_final
% % first assume stationary PM s.t. at each moment the forward generated
% % field is constant
% t_final = 10;  % [s]
% dt = 0.1;  % sampling time [s]
% t = 0:dt:t_final;  % time history
% n_t = length(t);  % length of time data
% 
% % Simulate measurement over time
% % Here consider circular array s.t. measurement is computed based on its
% % direction
% % B_meas_ideal rows: time index, columns: i-th sensor
% B_meas_ideal = ones(n_t,1)*B_measure_theoretical;
% 
% % Add noise into measurement on purpose with different level of magnitude
% noise_mu = 0;  % mean
% noise_sigma = 1e-7;  % variance
% noise = noise_mu + noise_sigma*randn(n_t,48);  % generating random noise according to given distribution
% 
% % Final corrupted measurement 
% B_meas = B_meas_ideal + noise;
% 
% % some plots
% figure
% plot(t,B_meas(:,6))
% hold on 
% plot([0,t(end)],[B_meas_ideal(1,6) B_meas_ideal(1,6)])
% xlabel('Time[s]')
% ylabel('Sensor measurement')
% legend('corrupted','ideal')

%% Backward estimation phase
% Given measurement of all sensors over a time horizon
% consider at one time, the measurement B_meas is 1 x n_sens vector

% Objective function: f_err = || B_est - B_meas ||^2
% Variable: pm_hat = (x_hat,y_hat,z_hat,theta_hat,phi_hat) as estimation of the position and orientation

% Disturbance field G
G_x = 0;
G_y = 0;
G_z = 0;
G = [G_x;G_y;G_z];

% initial guess on parameter estimation at time 0
pm_hat_0 = [-0.0;0.089;-0.00;pi/2;pi/2;M];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0)

% Compute Analytic Jacobian matrix： J <- R^(m x n) here
% J = Analytic_Jacobian_Mag_cir_array(pm_hat_0,xyz_s,meas_dir_s);


% % Levenberg–Marquardt least square algorithm solving for optimal position/orientation estimation
% % test
% B_meas = [0.002251179	-0.020839231	-0.171442474...
% 0.006875	0.142812392	-0.144833417...
% 0.01285366	0.161834937	0.105902026...
% 0.010088173	0.004057622	0.150368801...
% -0.015909321	0.004523515	-0.137550776...
% -0.014102385	0.145051269	-0.080955077...
% -0.002149897	0.125658074	0.136420404...
% 0.001431494	-0.011081308	0.152168353...
% 0.007019372	-0.003387205	-0.158306083...
% 0.117458167	0.241229167	-0.126916583...
% 0.057395353	0.119627699	0.161721744...
% 0.006132513	-0.01141294	0.134657994...
% -0.009514923	0.011241442	-0.133309808...
% -0.008916538	0.143974346	-0.054775538...
% -0.00439925	0.109192887	0.145673019...
% -0.00723566	-0.018559372	0.141078103...
% ]*1e-4;
% 
% % B_meas = [-0.054884641	-0.188673179	-1.630955385...
% % -0.035777833	1.444389167	-1.519999667...
% % -0.004972333	1.734195555	1.134027583...
% % 0.0009165	0.053888942	1.522195...
% % -0.163301314	0.072401708	-1.358174128...
% % -0.141640923	1.414640769	-0.891794769...
% % 0.039333462	1.299871015	1.267307308...
% % 0.079897462	-0.077589638	1.415100769...
% % 0.07445309	0.014440013	-1.390180827...
% % 0.303611083	1.373999833	-1.02980625...
% % 0.150749917	1.282334417	0.984083417...
% % 0.032611083	0.055166667	1.221527583...
% % -0.014692385	0.180356946	-1.154252115...
% % 0.019077077	1.2633323	-0.579076946...
% % 0.076948615	1.097820531	1.115537615...
% % 0.061307923	-0.0332051	1.214179...
% % ]*1e-4;


pm_hat_opt = PM_backward_estimation(B_meas',pm_hat_0,G,xyz_s,meas_dir_s); 

% weighting among parameters to be added!


%% evaluate estimated dipole model
% Disturbance field G
G_x = 0;
G_y = 0;
G_z = 0;
G = [G_x;G_y;G_z];

B_measure_theoretical = PM_forward_field(G,pm_hat_opt,xyz_s,meas_dir_s);
B_err = (B_meas - B_measure_theoretical)./B_measure_theoretical;
XYZ_estimation_plots(B_meas, B_measure_theoretical)
%% Timeline estimation
% use the current estimation for next step initial guess
% record the estimation history
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

% post data process: compute the standard deviation
mean_out = mean(pm_est_history);
sigma_out = std(pm_est_history);

% noise level ratio
J_1 = sigma_out/noise_sigma;

% position error


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







