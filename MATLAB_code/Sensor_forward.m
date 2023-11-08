function [B_meas] = Sensor_forward(pm_all,xyz_s,meas_dir_s,M)
% Read measurement of one sensor 3 directions for all pm posotions
% input:
% xyz_s: sensor array configuration <- 3x1
% pm_all: all PM positions&orientation <- 5xN
% meas_dir_s: normalized sensor measurement direction information <- R^(3 x 3)
% M: magnetic moment assumed

% output:
% B_meas: measured magnetic field depending on measurement direction (3*N x 1)
% B_sx: generated x component magnetic field
% B_sy: generated y component magnetic field
% B_sz: generated z component magnetic field

% extract parameters
[~,N] = size(pm_all);
[~,N_d] = size(meas_dir_s);
mu_0 = 4*pi*1e-7;
xyz_pm = pm_all(1:3,:);
theta0 = pm_all(4,:);
phi0 = pm_all(5,:);


% position vectors i-th column is i-th pm position vector
p_s =  xyz_s - xyz_pm; 
p_x_s = p_s(1,:);
p_y_s = p_s(2,:);
p_z_s = p_s(3,:);
p_s_sqr = diag(p_s'*p_s)';  % i-th column is x_i^2+y_i^2+z_i^2


% Forward magnetic field computations at all sensor locations
B_sx = mu_0/4/pi*M.*( 3*p_x_s.*(sin(theta0).*cos(phi0).*p_x_s + sin(theta0).*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - sin(theta0).*cos(phi0)./(p_s_sqr).^1.5 );
B_sy = mu_0/4/pi*M.*( 3*p_y_s.*(sin(theta0).*cos(phi0).*p_x_s + sin(theta0).*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - sin(theta0).*sin(phi0)./(p_s_sqr).^1.5 );
B_sz = mu_0/4/pi*M.*( 3*p_z_s.*(sin(theta0).*cos(phi0).*p_x_s + sin(theta0).*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - cos(theta0)./(p_s_sqr).^1.5 );

% project onto measurement direction
B_meas = meas_dir_s'*[B_sx;B_sy;B_sz];
B_meas = reshape(B_meas,[N_d*N,1]);


end