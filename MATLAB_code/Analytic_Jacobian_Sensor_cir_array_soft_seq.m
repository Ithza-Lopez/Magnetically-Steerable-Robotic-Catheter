%% Calibration with sensor soft iron distortion

%Experiment initialization
% choose magent type
clear;close all;clc;
col_labels = ["A","B","C"];
exp_date = '2-22-23';
magnet_type = 'A';
switch magnet_type
    case 'A'
        OD_mag = 3.175e-3;  % outer diameter of PM  [m]
        ID_mag = OD_mag/2;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.48; % magnetic remanance
        num_rows = -8:1:8;
    case 'B'
        V_pm = (1/4*25.4*1e-3)^3;
        Br = 1.48;
        num_rows = -8:1:8;
    case 'C'
        OD_mag = 3/16*25.4*1e-3;  % outer diameter of PM  [m]
        ID_mag = 1/16*25.4*1e-3;  % inner diameter of PM  [m]
        length_mag = OD_mag*2;   % length of PM  [m]
        V_pm = pi*(OD_mag^2-ID_mag^2^2)/4*length_mag;  % volume of PM
        Br = 1.32; % magnetic remanance
        num_rows = -6:1:6;
end

% choose row


% 4x4 Circular sensor array confuguration

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


%% individual evaluation for ith sensor

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

% M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
M_exp= 0.0396;

% chooose the number of sensor to calibrate as i
sensor = 1;

tilt_angle = 0; % 0 or 30

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,sensor);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(sensor-1)+1):3*sensor);

% measurement extraction
% choose the PM positions for calibration  
pm_all_all = [];
B_meas = [];
for i = 1:length(col_labels)
    col_letter = col_labels(i);
    for j = num_rows
        theta_pm = pi/2 ;  % orientation
        phi_pm = -pi/2- deg2rad(tilt_angle);  % orientation
        % extract pm positions
        pm_all = [Single_PM_Location(magnet_type, col_letter, j);theta_pm;phi_pm];
        pm_all_all = [pm_all_all pm_all];

        % collect measurement
        if tilt_angle >0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(j),'-', num2str(tilt_angle),'-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(j),'-Unbias.xlsx');
        end
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        B_j = Data((3*(sensor-1)+1):3*sensor)*1e-4;
        B_meas = [B_meas;B_j'];
    end
end

M_0 = M_exp;
[xyz_s_opt,M_opt] = Sensor_backward_estimation_noM(B_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_0);


N = length(B_meas)/3;

B_meas_sim = Sensor_forward(pm_all_all,xyz_s_opt,meas_dir_i_s_0,M_0);
B_sim_left = reshape(B_meas_sim',[N,3]);
B_meas_right = reshape(B_meas',[N,3]);

O = B_meas_right\B_sim_left;
O = O';




