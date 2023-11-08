function [J] = Analytic_Jacobian_Mag_sqr_array(M,pm0,xyz_s)
% This function computes the analytic Jacobian matrix only for square
% sensor array configuration, which means only Z component is utilized and
% computed 

% input:
% M: magnetic moment
% pm0: parameter in current step (x0,y0,z0,theta0,phi0) <- R^n  (more paremters later)
% xyz_s : sensor array information <- R^(3 x m)

% output:
% J: Jacobian matrix <- R^(m x n)

% General parameters
mu_0 = 4*pi*1e-7;
[~,m] = size(xyz_s);  % number of sensor
n = length(pm0);  % state dimension
p = xyz_s - pm0(1:3); % position vector for all sensors
theta0 = pm0(4);  % current theta
phi0 = pm0(5);  % current phi

% compute each element in Jacobian matrix
J = zeros(m,n);

% iterate for each sensor
for i = 1:m
    p_i = p(:,i);
    x_i = p_i(1);
    y_i = p_i(2);
    z_i = p_i(3);
    
    dB_iz_dx = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*z_i*( 4*x_i^2-y_i^2-z_i^2 )*sin(theta0)*cos(phi0) + ...
        15*x_i*y_i*z_i*sin(theta0)*sin(phi0) + 3*x_i*(-x_i^2-y_i^2+4*z_i^2)*cos(theta0) );
    
    dB_iz_dy = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 15*x_i*y_i*z_i*sin(theta0)*cos(phi0) + ...
        3*z_i*(-x_i^2+4*y_i^2-z_i^2)*sin(theta0)*sin(phi0) + 3*y_i*(-x_i^2-y_i^2+4*z_i^2)*cos(theta0) );
    
    dB_iz_dz = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*x_i*(-x_i^2-y_i^2+4*z_i^2)*sin(theta0)*cos(phi0) + ...
        3*y_i*(-x_i^2-y_i^2+4*z_i^2)*sin(theta0)*sin(phi0) + 3*z_i*(-3*x_i^2-3*y_i^2+2*z_i^2)*cos(theta0) );
    
    dB_iz_dtheta = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*( 3*x_i*z_i*cos(theta0)*cos(phi0) + ...
        3*y_i*z_i*cos(theta0)*sin(phi0) + (x_i^2+y_i^2-2*z_i^2)*sin(theta0) );
    
    dB_iz_dphi = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*(3*y_i*z_i*sin(theta0)*cos(phi0) - 3*x_i*z_i*sin(theta0)*sin(phi0) );
    
    J(i,:) = [ dB_iz_dx  dB_iz_dy  dB_iz_dz  dB_iz_dtheta  dB_iz_dphi];
end


end