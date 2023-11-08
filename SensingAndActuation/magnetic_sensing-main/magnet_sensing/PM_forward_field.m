function [B_meas] = PM_forward_field(M,G,pm_all,xyz_s,meas_dir_s,sat)
% input:
% M: magnetic moment value based on PM property (scalar)
% G: disturbance field 3x1
% xyz_s: sensor array configuration 3x16
% xyz_pm: current PM position 3x1
% meas_dir_s: normalized sensor measurement direction information <- R^(3 x m)
% sat: saturation magnetic field [T]

% output:
% B_meas: measured magnetic field depending on measurement direction (1 x m)
% B_sx: generated x component magnetic field
% B_sy: generated y component magnetic field
% B_sz: generated z component magnetic field

% extract parameters
mu_0 = 4*pi*1e-7;
xyz_pm = pm_all(1:3);
theta0 = pm_all(4);
phi0 = pm_all(5);

% Disturbance field G
G_x = G(1);
G_y = G(2);
G_z = G(3);

% position vectors i-th column is i-th sensor position vector
p_s = xyz_s - xyz_pm;
p_x_s = p_s(1,:);
p_y_s = p_s(2,:);
p_z_s = p_s(3,:);
p_s_sqr = diag(p_s'*p_s)';  % i-th column is x_i^2+y_i^2+z_i^2


% Forward magnetic field computations at all sensor locations
B_sx = G_x + mu_0/4/pi*M*( 3*p_x_s.*(sin(theta0)*cos(phi0).*p_x_s + sin(theta0)*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - sin(theta0)*cos(phi0)./(p_s_sqr).^1.5 );
B_sy = G_y + mu_0/4/pi*M*( 3*p_y_s.*(sin(theta0)*cos(phi0).*p_x_s + sin(theta0)*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - sin(theta0)*sin(phi0)./(p_s_sqr).^1.5 );
B_sz = G_z + mu_0/4/pi*M*( 3*p_z_s.*(sin(theta0)*cos(phi0).*p_x_s + sin(theta0)*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - cos(theta0)./(p_s_sqr).^1.5 );


% project onto measurement direction
B_meas = diag(meas_dir_s'*[B_sx;B_sy;B_sz])';

% consider saturation
if nargin >= 6
    B_meas = sign(B_meas).*min(abs(B_meas),sat);
end

end