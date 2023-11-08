function [B_meas] = PM_forward_field(pm_all,xyz_s,meas_dir_s)

% extract parameters
mu_0 = 4*pi*1e-7;
xyz_pm = pm_all(1:3);
theta0 = pm_all(4);
phi0 = pm_all(5);
M = pm_all(6);
[aa,m]=size(xyz_s);

p_s = xyz_s - repmat(xyz_pm,1,m);
p_x_s = p_s(1,:);
p_y_s = p_s(2,:);
p_z_s = p_s(3,:);
p_s_sqr = diag(p_s'*p_s)';  


B_sx = mu_0/4/pi*M*( 3*p_x_s.*(sin(theta0)*cos(phi0).*p_x_s + sin(theta0)*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - sin(theta0)*cos(phi0)./(p_s_sqr).^1.5 );
B_sy = mu_0/4/pi*M*( 3*p_y_s.*(sin(theta0)*cos(phi0).*p_x_s + sin(theta0)*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - sin(theta0)*sin(phi0)./(p_s_sqr).^1.5 );
B_sz = mu_0/4/pi*M*( 3*p_z_s.*(sin(theta0)*cos(phi0).*p_x_s + sin(theta0)*sin(phi0).*p_y_s + cos(theta0).*p_z_s)./(p_s_sqr).^2.5 - cos(theta0)./(p_s_sqr).^1.5 );


% project onto measurement direction
B_meas = diag(meas_dir_s'*[B_sx;B_sy;B_sz])';

end