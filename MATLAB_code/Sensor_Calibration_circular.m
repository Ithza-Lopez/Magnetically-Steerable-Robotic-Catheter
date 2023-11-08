%% Magnetically steerable catheter magnetic localization simulation
% model: ideal dipole model
% sensor configuration: 1. rectangular sensor array
%                       2. circular sensor array
% variables: distance 

% Start with 4x4 arrays in circular sensor plate array
clear;close all;clc;

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

load("B_measure.mat")
load('Sensor_of_interest.mat')

%% 4x4 Circular sensor array confuguration
clear;close all;clc;

n_sens = 16;

magnet_distance = 54.00; %mm
dim_repeat = magnet_distance; 

%all sensors working
% x_all = [repelem((dim_repeat),4), repelem((zeros(1,3)),4), repelem((-dim_repeat),4),  repelem((zeros(1,3)),4)]*1e-3; % 6 zeros and alternating x pos and neg
% y_all = [repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3),repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3),repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3), repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3)]*1e-3; % decreasing from 126.4 to 51.4 and increaseing to 126.4.
% z_all = [repelem((zeros(1,3)),4), repelem((dim_repeat),4),repelem((zeros(1,3)),4), repelem((-dim_repeat),4)]*1e-3; %alternating z pos and neg, and zeros
% xyz_all = [x_all;y_all;z_all];  % i-th column is i-th sensor location (x_i,y_i,z_i)

% index for all sensors ith column->ith sensor position
x_all = [repelem((dim_repeat),4), repelem((zeros(1,1)),4), repelem((-dim_repeat),4),  repelem((zeros(1,1)),4)]*1e-3; % 6 zeros and alternating x pos and neg
y_all = [repelem(-37.5,1),repelem(-12.5,1),repelem(12.5,1),repelem(37.5,1),repelem(-37.5,1),repelem(-12.5,1),repelem(12.5,1),repelem(37.5,1),repelem(-37.5,1),repelem(-12.5,1),repelem(12.5,1),repelem(37.5,1), repelem(-37.5,1),repelem(-12.5,1),repelem(12.5,1),repelem(37.5,1)]*1e-3; % decreasing from 126.4 to 51.4 and increaseing to 126.4.
z_all = [repelem((zeros(1,1)),4), repelem((dim_repeat),4),repelem((zeros(1,1)),4), repelem((-dim_repeat),4)]*1e-3; %alternating z pos and neg, and zeros
xyz_all = [x_all;y_all;z_all];



% define measurement direction vector (3 x m)
% i-th column shows the [x,y,z] components of the measurement at i-th
% sensor
% ith sensor's direction is 3(i-1)+1 ~ 3i
meas_dir_all_s= [fliplr(eye(3)).*[1;-1;1],...
              fliplr(eye(3)).*[1;-1;1],...
              fliplr(eye(3)).*[1;-1;1],...
              fliplr(eye(3)).*[1;-1;1],...
              eye(3).*[-1;-1;1],...
              eye(3).*[-1;-1;1],...
              eye(3).*[-1;-1;1],...
              eye(3).*[-1;-1;1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              fliplr(eye(3)).*[-1;-1;-1],...
              eye(3).*[1;-1;-1],...
              eye(3).*[1;-1;-1],...
              eye(3).*[1;-1;-1],...
              eye(3).*[1;-1;-1]];

meas_dir_all_s = meas_dir_all_s./sqrt(sum(meas_dir_all_s.^2));  % normalization to unit norm in each direction

%%  Calibration for all sensors 

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

% choose magent type
magnet_type = 'A';
switch magnet_type
    case 'A'
        OD_mag = 3.175e-3;  % outer diameter of PM  [m]
        ID_mag = OD_mag/2;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.48; % magnetic remanance
    case 'B'
        OD_mag = 3/16*25.4*1e-3;  % outer diameter of PM  [m]
        ID_mag = 1/16*25.4*1e-3;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.32; % magnetic remanance
    case 'C'
        V_pm = (1/4*25.4*1e-3)^3;
        Br = 1.48;
end

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
% M_exp = 0.0055;

% choose row
col_letter = 'B';

% chooose the number of sensor to calibrate as i
xyz_s_calibrated = [];
xyz_s_calibrated_all = [];
M_all = [];
for i = 1:16  % iterated among all sensors

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,i);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(i-1)+1):3*i);

% measurement extraction
% choose the PM positions for calibration
No_meas = -8:1:8;   
pm_all_all = [];
y_meas = [];

for j = No_meas
    theta_pm = pi/2;  % orientation
    phi_pm = -pi/2;  % orientation
    % extract pm positions
    pm_all = [Single_PM_Location(magnet_type, col_letter, j);theta_pm;phi_pm];
    pm_all_all = [pm_all_all pm_all]; 
    
    % collect measurement
    filename = ['2-9-23-A-B-',num2str(j),'-Unbias.xlsx'];
    Data = readmatrix(filename,'Sheet', 'Avg_vector');
    y_j = Data((3*(i-1)+1):3*i)*1e-4;
    y_meas = [y_meas;y_j'];
end 

% calibration
M_0 = M_exp;
[xyz_s_opt,M_opt] = Sensor_backward_estimation(y_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_0);

xyz_s_calibrated = [xyz_s_calibrated xyz_s_opt];
xyz_s_calibrated_all = [xyz_s_calibrated_all xyz_s_opt xyz_s_opt xyz_s_opt];
M_all = [M_all M_opt];

end



%% individual evaluation for ith sensor

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

% choose magent type
magnet_type = 'A';
switch magnet_type
    case 'A'
        OD_mag = 3.175e-3;  % outer diameter of PM  [m]
        ID_mag = OD_mag/2;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.48; % magnetic remanance
    case 'B'
        OD_mag = 3/16*25.4*1e-3;  % outer diameter of PM  [m]
        ID_mag = 1/16*25.4*1e-3;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.32; % magnetic remanance
    case 'C'
        V_pm = (1/4*25.4*1e-3)^3;
        Br = 1.48;
end

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment

% choose row
col_letter = 'B';

% chooose the number of sensor to calibrate as i
i = 1;

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,i);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(i-1)+1):3*i);

% measurement extraction
% choose the PM positions for calibration
No_meas = -8:1:8;   
pm_all_all = [];
y_meas = [];

for j = No_meas
    theta_pm = pi/2;  % orientation
    phi_pm = -pi/2;  % orientation
    % extract pm positions
    pm_all = [Single_PM_Location(magnet_type, col_letter, j);theta_pm;phi_pm];
    pm_all_all = [pm_all_all pm_all]; 
    
    % collect measurement
    filename = ['2-9-23-A-B-',num2str(j),'-Unbias.xlsx'];
    Data = readmatrix(filename,'Sheet', 'Avg_vector');
    y_j = Data((3*(i-1)+1):3*i)*1e-4;
    y_meas = [y_meas;y_j'];
end 

% calibration
M_0 = M_exp;
[xyz_s_opt,M_opt] = Sensor_backward_estimation(y_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_0);

y_measure_theoretical = Sensor_forward(pm_all_all,xyz_s_opt,meas_dir_i_s_0,M_opt);
y_err = (y_meas - y_measure_theoretical)./y_measure_theoretical;




%%  Redo estimation using new sensor layout
% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

% choose magent type
magnet_type = 'A';
switch magnet_type
    case 'A'
        OD_mag = 3.175e-3;  % outer diameter of PM  [m]
        ID_mag = OD_mag/2;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.48; % magnetic remanance
    case 'B'
        OD_mag = 3/16*25.4*1e-3;  % outer diameter of PM  [m]
        ID_mag = 1/16*25.4*1e-3;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.32; % magnetic remanance
    case 'C'
        V_pm = (1/4*25.4*1e-3)^3;
        Br = 1.48;
end

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment

% choose row
col_letter = 'B';

% choose the magnet to test
% No_mag = 1;
No_mag_all = -8:1:8;

pm_hat_0 = [-0.0;0.080;-0.00;pi/2;-pi/2;M_exp];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

pm_est_all = [];
pm_true_all = [];

for No_mag = No_mag_all
    
% PM orientation used in the test
pm_ori_true = [pi/2;-pi/2];

% extract ground truth position
pm_pos_true = Single_PM_Location(magnet_type, col_letter, No_mag);

% true parameter vector
pm_true = [pm_pos_true;pm_ori_true;M_exp];

pm_true_all = [pm_true_all pm_true];

% extract measurement for this specific pm position for all sensors
filename = ['2-9-23-A-B-',num2str(No_mag),'-Unbias.xlsx'];
Data = readmatrix(filename,'Sheet', 'Avg_vector');
B_meas = Data'*1e-4;  % B field measurement at all sensors

pm_hat_opt = PM_backward_estimation(B_meas,pm_hat_0,xyz_s_calibrated_all,meas_dir_all_s); 

pm_est_all = [pm_est_all pm_hat_opt];

% use previous true position as a initial guess
pm_hat_0 = pm_true;   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

end


%% 
figure(1)
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend('True','Estimation')
xlim([-20 2])
ylim([-90 90])

figure(2)
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend('True','Estimation')
xlim([-2 2])
ylim([-90 90])


