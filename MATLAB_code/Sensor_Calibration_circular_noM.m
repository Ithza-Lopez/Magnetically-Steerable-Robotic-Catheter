%% Magnetically steerable catheter magnetic localization simulation
%Experiment initialization
% choose magent type
clear;close all;clc;
col_labels = ["B","C", "D"];
exp_date = '5-16-23';
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

%%  Calibration for all sensors using only straight data or tilted data

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;
% M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
tilt_angle = 0; % 0 or 30
% M_exp = 0.052;
% M_exp = 0.0396;
% M_exp = 0.0358;]
M_exp = 0.0395266;

num_rows = -8:1:8;

% chooose the number of sensor to calibrate as i
xyz_s_calibrated = [];
xyz_s_calibrated_all = [];
M_all = [];
for sensor = 1:n_sens  % iterated among all sensors

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,sensor);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(sensor-1)+1):3*sensor);

% measurement extraction
% choose the PM positions for calibration
% height_range = [-20,0,20];
height_range = [0];
pm_all_all = [];
y_meas = [];
for height = 1:length(height_range)
    height_num = height_range(height);
for col = 1:length(col_labels)
    col_letter = col_labels(col);
    for row = num_rows
        theta_pm = pi/2  ;  % orientation
        phi_pm = -pi/2- deg2rad(tilt_angle);  % orientation
        % extract pm positions
        pm_all = [Single_PM_Location(magnet_type, col_letter, row,0);theta_pm;phi_pm];
        pm_all_all = [pm_all_all pm_all];

        % collect measurement
        filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-',num2str(height_num), '-Unbias.xlsx');
%         folder_name = 'Middle_Plate';
%         file_path = fullfile(pwd, folder_name, filename);
%         Data = readmatrix(file_path,'Sheet', 'Avg_vector');
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        y_j = Data((3*(sensor-1)+1):3*sensor)*1e-4;
        y_meas = [y_meas;y_j'];
    end
end
end 
% calibration
M_0 = M_exp;
[xyz_s_opt,M_opt] = Sensor_backward_estimation_noM(y_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_0);

xyz_s_calibrated = [xyz_s_calibrated xyz_s_opt];
xyz_s_calibrated_all = [xyz_s_calibrated_all xyz_s_opt xyz_s_opt xyz_s_opt];
M_all = [M_all M_opt];

end

disp('Completed section')


%%  Calibration for all sensors using all data

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;
% M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
% tilt_angle = 30;
% M_exp = 0.0396;
M_exp = 0.0395266;

num_rows = -4:1:4;

% chooose the number of sensor to calibrate as i
xyz_s_calibrated = [];
xyz_s_calibrated_all = [];
M_all = [];

tilt_angle_all = [0;30];



for sensor = 1:n_sens  % iterated among all sensors

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,sensor);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(sensor-1)+1):3*sensor);

% measurement extraction
% choose the PM positions for calibration
   
pm_all_all = [];
y_meas = [];
for iii = 1:length(tilt_angle_all)
    tilt_angle = tilt_angle_all(iii);
for col = 1:length(col_labels)
    col_letter = col_labels(col);
    for row = num_rows
        theta_pm = pi/2  ;  % orientation
        phi_pm = -pi/2- deg2rad(tilt_angle);  % orientation
        % extract pm positions
        pm_all = [Single_PM_Location(magnet_type, col_letter, row);theta_pm;phi_pm];
        pm_all_all = [pm_all_all pm_all];

        % collect measurement
        if tilt_angle ~=0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-Unbias.xlsx');
        end
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        y_j = Data((3*(sensor-1)+1):3*sensor)*1e-4;
        y_meas = [y_meas;y_j'];
    end
end
end

% calibration
M_0 = M_exp;
[xyz_s_opt,M_opt] = Sensor_backward_estimation_noM(y_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_0);

xyz_s_calibrated = [xyz_s_calibrated xyz_s_opt];
xyz_s_calibrated_all = [xyz_s_calibrated_all xyz_s_opt xyz_s_opt xyz_s_opt];
M_all = [M_all M_opt];

end




%% individual evaluation for ith sensor

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment

% chooose the number of sensor to calibrate as i
sensor = 1;

tilt_angle = 0; % 0 or 30

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,sensor);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(sensor-1)+1):3*sensor);

% measurement extraction
% choose the PM positions for calibration  
pm_all_all = [];
y_meas = [];
for col = 1:length(col_labels)
    col_letter = col_labels(col);
    for row = num_rows
        theta_pm = pi/2 ;  % orientation
        phi_pm = -pi/2- deg2rad(tilt_angle);  % orientation
        % extract pm positions
        pm_all = [Single_PM_Location(magnet_type, col_letter, row,0);theta_pm;phi_pm];
        pm_all_all = [pm_all_all pm_all];

        % collect measurement
        if tilt_angle >0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-Unbias.xlsx');
        end
        folder_name = 'Middle_Plate';
        file_path = fullfile(pwd, folder_name, filename);
        Data = readmatrix(file_path,'Sheet', 'Avg_vector');
        y_j = Data((3*(sensor-1)+1):3*sensor)*1e-4;
        y_meas = [y_meas;y_j'];
    end
end
% calibration
M_0 = M_exp;
[xyz_s_opt,M_opt] = Sensor_backward_estimation_noM(y_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_0);

y_measure_theoretical = Sensor_forward(pm_all_all,xyz_s_opt,meas_dir_i_s_0,M_opt);
y_err = (y_meas - y_measure_theoretical)./y_measure_theoretical;




%%  Redo estimation using new sensor layout
% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

% M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
% M_exp = 0.0396;
% M_exp = 0.0358;
M_exp = 0.0395266;
% choose the magnet to test

% pm_hat_0 = [-0.0;0.080;-0.00;pi/2;-pi/2- deg2rad(tilt_angle);M_exp];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

pm_est_all = [];
pm_true_all = [];
tilt_angle = 30;
num_rows = -8:1:8;
for col = 1:length(col_labels)
    col_letter = col_labels(col);
    for row = num_rows

        % PM orientation used in the test
        pm_ori_true = [pi/2;-pi/2-deg2rad(tilt_angle)];

        % extract ground truth position
        pm_pos_true = Single_PM_Location(magnet_type, col_letter, row,0);

        % true parameter vector
%         pm_true = [pm_pos_true;pm_ori_true;M_exp];
        pm_true = [pm_pos_true;pm_ori_true];
        pm_hat_0 = pm_true+[2e-3;2e-3;2e-3;0.01;-0.01];

        pm_true_all = [pm_true_all pm_true];

        % extract measurement for this specific pm position for all sensors
        if tilt_angle ~=0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-0-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-Unbias.xlsx');
        end
%         folder_name = 'Middle_Plate';
%         file_path = fullfile(pwd, folder_name, filename);
%         Data = readmatrix(file_path,'Sheet', 'Avg_vector');
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        B_meas = Data'*1e-4;  % B field measurement at all sensors

        % use true position as a initial guess
%         pm_hat_0 = pm_true;   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)
        pm_hat_opt = PM_backward_estimation_noM(B_meas,pm_hat_0,xyz_s_calibrated_all,meas_dir_all_s,M_exp);

        pm_est_all = [pm_est_all pm_hat_opt];



    end
end
% compute error in [mm]
pos_error = (pm_est_all(1:3,:)-pm_true_all(1:3,:))*1e3;
ang_error = (pm_est_all(4:5,:)-pm_true_all(4:5,:))/pi*180;
pos_sqr_err = sqrt(sum(pos_error.^2,1));
N_D = length(pos_sqr_err);
pos_sum_sqr_err = sum(pos_sqr_err)/N_D;
theta_rms_error = sum(abs(ang_error(1,:)))/N_D;
phi_rms_error = sum(abs(ang_error(2,:)))/N_D;

save('Angled_Cal_pm_est_all.mat','pm_est_all' )
save('pm_true_all.mat','pm_true_all')

%%  Creating magnet estimation figures
close all

% figure(1)
% plot3(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
% hold on
% plot3(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3, pm_est_all(3,:)*1e3,'rx')
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('z[mm]')
% legend('True','Estimation')
% % xlim([-2 2])
% % ylim([-90 90])
% % zlim([-2 2])

figure(2)
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend('True','Estimation')
% xlim([-2 2])
% ylim([-90 90])


figure(3)
plot(pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
hold on
plot(pm_est_all(2,:)*1e3,pm_est_all(3,:)*1e3,'rx')
xlabel('y[mm]')
ylabel('z[mm]')
legend('True','Estimation')
% xlim([-90 90])
% ylim([-2.5 2.5])

% saveas(figure(1), 'AfterCal XYZ View.jpeg' )
% saveas(figure(2), 'AfterCal Top View.jpeg' )
% saveas(figure(3), 'AfterCal Side View.jpeg' )
% XYZ_estimation_plots(pm_true_all, pm_est_all)
% saveas(figure(4), 'AfterCal XYZ_estimation_plots.jpeg' )
%% evaluate new estimated dipole model

B_measure_theoretical = PM_forward_field(pm_hat_opt,xyz_s_calibrated_all,meas_dir_all_s);
B_err = abs(B_meas' - B_measure_theoretical);
B_err_ratio = B_err./B_meas';




%% Grinding M

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;
M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
% tilt_angle = 30;
% M_exp = 0.0055;

num_rows = -8:1:8;

% chooose the number of sensor to calibrate as i
xyz_s_calibrated = [];
xyz_s_calibrated_all = [];
% M_all = [];

tilt_angle_all = [0];

M_all = M_exp*6/10:5e-4:M_exp*12/10;
y_err_all = [];
y_err_norm_all = [];

for k = 1:length(M_all)
k
M_k = M_all(k);
y_err = 0;



for sensor = 1:n_sens  % iterated among all sensors

% extract its position and orientation  as initial guess
xyz_i_s_0 = xyz_all(:,sensor);
meas_dir_i_s_0 = meas_dir_all_s(:,(3*(sensor-1)+1):3*sensor);

% measurement extraction
% choose the PM positions for calibration
   
pm_all_all = [];
y_meas = [];
for iii = 1:length(tilt_angle_all)
    tilt_angle = tilt_angle_all(iii);
for col = 1:length(col_labels)
    col_letter = col_labels(col);
    for row = num_rows
        theta_pm = pi/2  ;  % orientation
        phi_pm = -pi/2- deg2rad(tilt_angle);  % orientation
        % extract pm positions
        pm_all = [Single_PM_Location(magnet_type, col_letter, row, 0);theta_pm;phi_pm];
        pm_all_all = [pm_all_all pm_all];

        % collect measurement
        if tilt_angle ~=0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-',num2str(tilt_angle),'-0','-Unbias.xlsx');
        end
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        y_j = Data((3*(sensor-1)+1):3*sensor)*1e-4;
        y_meas = [y_meas;y_j'];
    end 
end
end

% calibration
[xyz_s_opt,M_opt] = Sensor_backward_estimation_noM(y_meas,pm_all_all,xyz_i_s_0,meas_dir_i_s_0,M_k);

% xyz_s_calibrated = [xyz_s_calibrated xyz_s_opt];
% xyz_s_calibrated_all = [xyz_s_calibrated_all xyz_s_opt xyz_s_opt xyz_s_opt];
% M_all = [M_all M_opt];
y_measure_theoretical = Sensor_forward(pm_all_all,xyz_s_opt,meas_dir_i_s_0,M_k);
y_err = y_err + abs(y_meas - y_measure_theoretical);


end

y_err_all = [y_err_all y_err];
y_err_norm_all = [y_err_norm_all norm(y_err,2)];

end

%% plots
figure()
plot(M_all,y_err_norm_all)
xlabel('M')
ylabel('Error Norm')
min_error = min(y_err_norm_all);
M_cal = find(y_err_norm_all == min(y_err_norm_all))
%





