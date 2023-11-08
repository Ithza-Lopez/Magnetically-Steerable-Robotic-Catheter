%% Magnetically steerable catheter magnetic localization simulation
% model: ideal dipole model

% variables: distance 

% Start with 4x4 arrays in circular sensor plate array
clear;close all;clc;

exp_date = '5-16-23';
magnet_type = 'A';
col_labels = ["B","C","D"];
% col_labels = ["C"];

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;


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


% 4x4 Circular sensor array confuguration

n_sens = 16;

magnet_distance = 54; %mm
dim_repeat = ones(1,3)*magnet_distance; 

%all sensors working
% x_s = [repelem((zeros(1,3)),4), repelem((-dim_repeat),4), repelem((dim_repeat),4), repelem((zeros(1,3)),4)]*1e-3; % 6 zeros and alternating x pos and neg
% y_s = [repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3),repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3),repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3),repelem(51.4,3),repelem(76.4,3),repelem(101.4,3),repelem(126.4,3)]*1e-3-0.0889; % decreasing from 126.4 to 51.4 and increaseing to 126.4.
% z_s = [repelem((dim_repeat),4),repelem((zeros(1,3)),8), repelem((-dim_repeat),4)]*1e-3; %alternating z pos and neg, and zeros
% xyz_s = [x_s;y_s;z_s];  % i-th column is i-th sensor location (x_i,y_i,z_i)

x_all = [repelem((dim_repeat),4), repelem((zeros(1,3)),4), repelem((-dim_repeat),4),  repelem((zeros(1,3)),4)]*1e-3; % 6 zeros and alternating x pos and neg
y_all = [repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3),repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3),repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3), repelem(-37.5,3),repelem(-12.5,3),repelem(12.5,3),repelem(37.5,3)]*1e-3; % decreasing from 126.4 to 51.4 and increaseing to 126.4.
z_all = [repelem((zeros(1,3)),4), repelem((dim_repeat),4),repelem((zeros(1,3)),4), repelem((-dim_repeat),4)]*1e-3; %alternating z pos and neg, and zeros
xyz_s = [x_all;y_all;z_all];  % i-th column is i-th sensor location (x_i,y_i,z_i)



% define measurement direction vector (3 x m)
% i-th column shows the [x,y,z] components of the measurement at i-th
% sensor
%all sensor working
meas_dir_s= [fliplr(eye(3)).*[1;-1;1],...
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

meas_dir_s = meas_dir_s./sqrt(sum(meas_dir_s.^2));  % normalization to unit norm in each direction

%% Test all positions and report/plot all results
% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

% M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
M_exp = 0.0395
theta = pi/2;
phi = -pi/2;

%Change these 2 for new columns and heights
tilt_angle = 30;
pm_hat_0 = [-0.020;0.080;0.000;theta; phi-deg2rad(tilt_angle);M_exp];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

pm_est_all = [];
pm_true_all = [];

height_range = [0];
num_rows = -8:1:8;
for height = 1:length(height_range)
    height_num = height_range(height);
for i = 1:length(col_labels)
    col_letter = col_labels(i);
    for row = num_rows
        
        % PM orientation used in the test
        pm_ori_true = [theta;phi-deg2rad(tilt_angle)];

        % extract ground truth position
        pm_pos_true = Single_PM_Location(magnet_type, col_letter, row, height_num);

        % true parameter vector
        pm_true = [pm_pos_true;pm_ori_true;M_exp];
           
        pm_true_all = [pm_true_all pm_true];

        % extract measurement for this specific pm position for all sensors
        filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-',num2str(height_num), '-Unbias.xlsx');
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        B_meas = Data'*1e-4;  % B field measurement at all sensors
        pm_hat_0 = pm_true;   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

        pm_hat_opt = PM_backward_estimation(B_meas,pm_hat_0,xyz_s,meas_dir_s);

        pm_est_all = [pm_est_all pm_hat_opt];
        % use previous true position as a initial guess
        
    end
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

save('No_cal_pm_est_all.mat','pm_est_all' )
% save('pm_true_all.mat','pm_true_all')
%% 
close all
% figure()
% plot3(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3, pm_true_all(3,:)*1e3,'bo')
% hold on
% plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3, pm_est_all(3,:)*1e3,'rx')
% xlabel('x[mm]')
% ylabel('y[mm]')
% legend('True','Estimation')
% xlim([-20 2])
% ylim([-90 90])

figure(1)
plot3(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
hold on
plot3(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3, pm_est_all(3,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('z[mm]')
legend('True','Estimation')
title('XYZ View')
% xlim([-2 2])
% ylim([-90 90])
% zlim([-2 2])


figure(2)
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend('True','Estimation')
title('Top View')
% xlim([-25 25])
% ylim([-50 50])
% axis equal

figure(3)
plot(pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
hold on
plot(pm_est_all(2,:)*1e3,pm_est_all(3,:)*1e3,'rx')
xlabel('y[mm]')
ylabel('z[mm]')
legend('True','Estimation')
title('Side View')
% % xlim([-90 90])
% % ylim([-2.5 2.5])

% saveas(figure(1), 'NoCal XYZ View.jpeg' )
% saveas(figure(2), 'NoCal Top View.jpeg' )
% saveas(figure(3), 'NoCal Side View.jpeg' )
% XYZ_estimation_plots(pm_true_all, pm_est_all)
% saveas(figure(4), 'NoCal XYZ_estimation_plots.jpeg' )


%% Test individual position
% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
theta = pi/2;
tilt_angle = 30;
phi = -pi/2-deg2rad(tilt_angle);


pm_hat_0 = [-0.00;0.00;-0.00;theta;phi;M_exp];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

col_letter = 'A';
row = 0;
        
% PM orientation used in the test
pm_ori_true = [theta;phi];

% extract ground truth position
pm_pos_true = Single_PM_Location(magnet_type, col_letter, row);

% true parameter vector
pm_true = [pm_pos_true;pm_ori_true;M_exp];

% extract measurement for this specific pm position for all sensors

if tilt_angle ~=0
    filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-Unbias.xlsx');
else
    filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-Unbias.xlsx');
end
Data = readmatrix(filename,'Sheet', 'Avg_vector');
B_meas = Data'*1e-4;  % B field measurement at all sensors
% pm_hat_0 = pm_true;   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

pm_hat_opt = PM_backward_estimation(B_meas,pm_hat_0,xyz_s,meas_dir_s);

%% evaluate estimated dipole model

B_measure_theoretical = PM_forward_field(pm_hat_opt,xyz_s,meas_dir_s);
B_measure_theoretical = B_measure_theoretical';
B_err = abs(B_meas - B_measure_theoretical);
B_err_ratio = B_err./B_meas;


%%
mu_0 = 4*pi*1e-7;

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
theta = pi/2;
phi = -pi/2;
tilt_angle = 30;

pm_hat_0 = [-0.010;0.080;-0.00;theta;-phi-deg2rad(tilt_angle);M_exp];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

pm_est_all = [];
pm_true_all = [];

num_rows = -8:1:8;

for i = 1:length(col_labels)
    col_letter = col_labels(i);
    for row = num_rows
        
        % PM orientation used in the test
        pm_ori_true = [theta;phi-deg2rad(tilt_angle)];

        % extract ground truth position
        pm_pos_true = Single_PM_Location(magnet_type, col_letter, row);

        % true parameter vector
        pm_true = [pm_pos_true;pm_ori_true;M_exp];
           
        pm_true_all = [pm_true_all pm_true];

        % extract measurement for this specific pm position for all sensors

        if tilt_angle ~=0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-Unbias.xlsx');
        end
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        B_meas = Data'*1e-4;  % B field measurement at all sensors
        pm_hat_0 = pm_true;   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)
        threshold = 1e-6;

        % weight inversely proportional to distance
        Dist = sum((xyz_s - pm_pos_true).^2,1);
        W = diag((1./Dist).*(abs(B_meas')>threshold));

        pm_hat_opt = PM_backward_estimation_WLS(B_meas,pm_hat_0,xyz_s,meas_dir_s,W); 

        pm_est_all = [pm_est_all pm_hat_opt];
        % use previous true position as a initial guess
        
    end
end

% compute error in [mm]
pos_error = (pm_est_all(1:3,:)-pm_true_all(1:3,:))*1e3;
ang_error = (pm_est_all(4:5,:)-pm_true_all(4:5,:))/pi*180;
pos_sqr_err = sqrt(sum(pos_error.^2,1));




%%


figure(1)
plot3(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
hold on
plot3(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3, pm_est_all(3,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('z[mm]')
legend('True','Estimation')
title('XYZ View')
% xlim([-2 2])
% ylim([-90 90])
% zlim([-2 2])


figure(2)
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend('True','Estimation')
title('Top View')
% xlim([-2 2])
% ylim([-90 90])


figure(3)
plot(pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
hold on
plot(pm_est_all(2,:)*1e3,pm_est_all(3,:)*1e3,'rx')
xlabel('y[mm]')
ylabel('z[mm]')
legend('True','Estimation')
title('Side View')



%%   Before calibration, localization without updating M

% Magnetic moment M computation
mu_0 = 4*pi*1e-7;

M_exp= V_pm*Br/mu_0;  % estimated magnetic moment
% choose the magnet to test

% pm_hat_0 = [-0.0;0.080;-0.00;pi/2;-pi/2- deg2rad(tilt_angle);M_exp];   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)

pm_est_all = [];
pm_true_all = [];
tilt_angle = 00;
num_rows = -4:1:4;
for i = 1:length(col_labels)
    col_letter = col_labels(i);
    for j = num_rows

        % PM orientation used in the test
        pm_ori_true = [pi/2;-pi/2-deg2rad(tilt_angle)];

        % extract ground truth position
        pm_pos_true = Single_PM_Location(magnet_type, col_letter, j);

        % true parameter vector
        pm_true = [pm_pos_true;pm_ori_true];
        pm_hat_0 = pm_true;

        pm_true_all = [pm_true_all pm_true];

        % extract measurement for this specific pm position for all sensors
        if tilt_angle ~=0
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(j),'-', num2str(tilt_angle),'-Unbias.xlsx');
        else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(j),'-Unbias.xlsx');
        end
        Data = readmatrix(filename,'Sheet', 'Avg_vector');
        B_meas = Data'*1e-4;  % B field measurement at all sensors

        % use true position as a initial guess
        pm_hat_0 = pm_true;   % (x_hat0,y_hat0,z_hat0,theta_hat0,phi_hat0,M0)
        pm_hat_opt = PM_backward_estimation_noM(B_meas,pm_hat_0,xyz_s,meas_dir_s,M_exp);

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
% 
% if tilt_angle == 0 
%     save('Straight_No_cal_pm_est_all.mat','pm_est_all' )
%     save('Straight_pm_true_all','pm_true_all')
% else
%     save('Angle_No_cal_pm_est_all.mat','pm_est_all' )
%     save('Angle_pm_true_all','pm_true_all')
% end

%%


% % figure(1)
% % plot3(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
% % hold on
% % plot3(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3, pm_est_all(3,:)*1e3,'rx')
% % xlabel('x[mm]')
% % ylabel('y[mm]')
% % zlabel('z[mm]')
% % legend('True','Estimation')
% % title('XYZ View')
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
title('Top View')
% xlim([-2 2])
% ylim([-90 90])


% figure(3)
% plot(pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
% hold on
% plot(pm_est_all(2,:)*1e3,pm_est_all(3,:)*1e3,'rx')
% xlabel('y[mm]')
% ylabel('z[mm]')
% legend('True','Estimation')
% title('Side View')

saveas(figure(2), 'Angle Uncalibrated Top View.pdf' )















