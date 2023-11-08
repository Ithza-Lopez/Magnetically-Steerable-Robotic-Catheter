function [J] = Analytic_Jacobian_Mag_cir_array(M,pm0,xyz_s,meas_dir_s)
% This function computes the analytic Jacobian matrix for circular tube
% sensor array configuration, which means the scalar measurement directions can
% be different but given by meas_dir_s matrix.

% input:
% M: magnetic moment
% pm0: parameter in current step (x0,y0,z0,theta0,phi0) <- R^n  (more paremters later)
% xyz_s : sensor array information <- R^(3 x m)
% meas_dir_s: normalized sensor measurement direction information <- R^(3 x m)

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
    % position components at i-th sensor
    p_i = p(:,i);
    x_i = p_i(1);
    y_i = p_i(2);
    z_i = p_i(3);
    
    % measurement direction at i-th sensor
    uvw_i = meas_dir_s(:,i);
    u_i = uvw_i(1);
    v_i = uvw_i(2);
    w_i = uvw_i(3);
    
    dB_ix_dx = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*x_i*( 2*x_i^2-3*y_i^2-3*z_i^2 )*sin(theta0)*cos(phi0) + ...
        3*y_i*( 4*x_i^2-y_i^2-z_i^2 )*sin(theta0)*sin(phi0) + 3*z_i*( 4*x_i^2-y_i^2-z_i^2 )*cos(theta0) );  
    
    dB_iy_dx = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*y_i*( 4*x_i^2-y_i^2-z_i^2 )*sin(theta0)*cos(phi0) + ...
        3*x_i*( -x_i^2+4*y_i^2-z_i^2 )*sin(theta0)*sin(phi0) + 15*x_i*y_i*z_i*cos(theta0) );
    
    dB_iz_dx = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*z_i*( 4*x_i^2-y_i^2-z_i^2 )*sin(theta0)*cos(phi0) + ...
        15*x_i*y_i*z_i*sin(theta0)*sin(phi0) + 3*x_i*( -x_i^2-y_i^2+4*z_i^2 )*cos(theta0) );
    
    dB_ix_dy = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*y_i*( 4*x_i^2-y_i^2-z_i^2 )*sin(theta0)*cos(phi0) + ...
        3*x_i*( -x_i^2+4*y_i^2-z_i^2 )*sin(theta0)*sin(phi0) + 15*x_i*y_i*z_i*cos(theta0) );
    
    dB_iy_dy = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*x_i*( -x_i^2+4*y_i^2-z_i^2 )*sin(theta0)*cos(phi0) + ...
        3*y_i*( -3*x_i^2+2*y_i^2-3*z_i^2 )*sin(theta0)*sin(phi0) + 3*z_i*( -x_i^2+4*y_i^2-z_i^2 )*cos(theta0)   );
    
    dB_iz_dy = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 15*x_i*y_i*z_i*sin(theta0)*cos(phi0) + ...
        3*z_i*( -x_i^2+4*y_i^2-z_i^2 )*sin(theta0)*sin(phi0) + 3*y_i*( -x_i^2-y_i^2+4*z_i^2 )*cos(theta0) );
    
    dB_ix_dz = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*z_i*( 4*x_i^2-y_i^2-z_i^2 )*sin(theta0)*cos(phi0) + ...
        15*x_i*y_i*z_i*sin(theta0)*sin(phi0) + 3*x_i*( -x_i^2-y_i^2+4*z_i^2 )*cos(theta0) );
    
    dB_iy_dz = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 15*x_i*y_i*z_i*sin(theta0)*cos(phi0) + ...
        3*z_i*( -x_i^2+4*y_i^2-z_i^2 )*sin(theta0)*sin(phi0) + 3*y_i*( -x_i^2-y_i^2+4*z_i^2 )*cos(theta0)  );
    
    dB_iz_dz = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-3.5*( 3*x_i*(-x_i^2-y_i^2+4*z_i^2)*sin(theta0)*cos(phi0) + ...
        3*y_i*( -x_i^2-y_i^2+4*z_i^2 )*sin(theta0)*sin(phi0) + 3*z_i*( -3*x_i^2-3*y_i^2+2*z_i^2 )*cos(theta0) );
    
    dB_ix_dtheta = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*( ( 2*x_i^2-y_i^2-z_i^2 )*cos(theta0)*cos(phi0) +...
        3*x_i*y_i*cos(theta0)*sin(phi0) -3*x_i*z_i*sin(theta0) );
    
    dB_iy_dtheta = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*( 3*x_i*y_i*cos(theta0)*cos(phi0) +...
        (-x_i^2+2*y_i^2-z_i^2)*cos(theta0)*sin(phi0) - 3*y_i*z_i*sin(theta0) );
    
    dB_iz_dtheta = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*( 3*x_i*z_i*cos(theta0)*cos(phi0) + ...
        3*y_i*z_i*cos(theta0)*sin(phi0) + ( x_i^2+y_i^2-2*z_i^2 )*sin(theta0) );
    
    dB_ix_dphi = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*( 3*x_i*y_i*sin(theta0)*cos(phi0) - ( 2*x_i^2-y_i^2-z_i^2)*sin(theta0)*sin(phi0) );
    
    dB_iy_dphi = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*( (-x_i^2+2*y_i^2-z_i^2)*sin(theta0)*cos(phi0) -3*x_i*y_i*sin(theta0)*sin(phi0)  );
    
    dB_iz_dphi = mu_0*M/4/pi*(x_i^2+y_i^2+z_i^2)^-2.5*(3*y_i*z_i*sin(theta0)*cos(phi0) - 3*x_i*z_i*sin(theta0)*sin(phi0) );
    
    dB_mi_dx = u_i*dB_ix_dx + v_i*dB_iy_dx + w_i*dB_iz_dx;
    
    dB_mi_dy = u_i*dB_ix_dy + v_i*dB_iy_dy + w_i*dB_iz_dy;
    
    dB_mi_dz = u_i*dB_ix_dz + v_i*dB_iy_dz + w_i*dB_iz_dz;
    
    dB_mi_dtheta = u_i*dB_ix_dtheta + v_i*dB_iy_dtheta + w_i*dB_iz_dtheta;
    
    dB_mi_dphi = u_i*dB_ix_dphi + v_i*dB_iy_dphi + w_i*dB_iz_dphi;
    
    J(i,:) = [ dB_mi_dx  dB_mi_dy  dB_mi_dz  dB_mi_dtheta  dB_mi_dphi];
end




end