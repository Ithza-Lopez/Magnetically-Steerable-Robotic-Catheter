function [J] = Analytic_Jacobian_Mag_cir_array_test(pm0,xyz_s,meas_dir_s)

mu_0 = 4*pi*1e-7;
[aaa,m] = size(xyz_s);  % number of sensor
p = xyz_s - repmat(pm0(1:3),1,m); % position vector for all sensors
theta0 = pm0(4);  % current thetaxu
phi0 = pm0(5);  % current phi
M = pm0(6);  % current magnetic moment

% position components at i-th sensor
x = p(1,:);
y = p(2,:);
z = p(3,:);

% pre-computations
xyz_sqr = sum(p.^2,1);
xyz_35 = xyz_sqr.^-3.5;
xyz_25 = xyz_sqr.^-2.5;

dB_ix_dx = mu_0*M/4/pi*xyz_35.*( 3*x.*( 2*x.^2-3*y.^2-3*z.^2 )*sin(theta0)*cos(phi0) + ...
    3*y.*( 4*x.^2-y.^2-z.^2 )*sin(theta0)*sin(phi0) + 3*z.*( 4*x.^2-y.^2-z.^2 )*cos(theta0) );

dB_iy_dx = mu_0*M/4/pi*xyz_35.*( 3*y.*( 4*x.^2-y.^2-z.^2 )*sin(theta0)*cos(phi0) + ...
    3*x.*( -x.^2+4*y.^2-z.^2 )*sin(theta0)*sin(phi0) + 15*x.*y.*z*cos(theta0) );

dB_iz_dx = mu_0*M/4/pi*xyz_35.*( 3*z.*( 4*x.^2-y.^2-z.^2 )*sin(theta0)*cos(phi0) + ...
    15*x.*y.*z*sin(theta0)*sin(phi0) + 3*x.*( -x.^2-y.^2+4*z.^2 )*cos(theta0) );

dB_ix_dy = mu_0*M/4/pi*xyz_35.*( 3*y.*( 4*x.^2-y.^2-z.^2 )*sin(theta0)*cos(phi0) + ...
    3*x.*( -x.^2+4*y.^2-z.^2 )*sin(theta0)*sin(phi0) + 15*x.*y.*z*cos(theta0) );

dB_iy_dy = mu_0*M/4/pi*xyz_35.*( 3*x.*( -x.^2+4*y.^2-z.^2 )*sin(theta0)*cos(phi0) + ...
    3*y.*( -3*x.^2+2*y.^2-3*z.^2 )*sin(theta0)*sin(phi0) + 3*z.*( -x.^2+4*y.^2-z.^2 )*cos(theta0)   );

dB_iz_dy = mu_0*M/4/pi*xyz_35.*( 15*x.*y.*z*sin(theta0)*cos(phi0) + ...
    3*z.*( -x.^2+4*y.^2-z.^2 )*sin(theta0)*sin(phi0) + 3*y.*( -x.^2-y.^2+4*z.^2 )*cos(theta0) );

dB_ix_dz = mu_0*M/4/pi*xyz_35.*( 3*z.*( 4*x.^2-y.^2-z.^2 )*sin(theta0)*cos(phi0) + ...
    15*x.*y.*z*sin(theta0)*sin(phi0) + 3*x.*( -x.^2-y.^2+4*z.^2 )*cos(theta0) );

dB_iy_dz = mu_0*M/4/pi*xyz_35.*( 15*x.*y.*z*sin(theta0)*cos(phi0) + ...
    3*z.*( -x.^2+4*y.^2-z.^2 )*sin(theta0)*sin(phi0) + 3*y.*( -x.^2-y.^2+4*z.^2 )*cos(theta0)  );

dB_iz_dz = mu_0*M/4/pi*xyz_35.*( 3*x.*(-x.^2-y.^2+4*z.^2)*sin(theta0)*cos(phi0) + ...
    3*y.*( -x.^2-y.^2+4*z.^2 )*sin(theta0)*sin(phi0) + 3*z.*( -3*x.^2-3*y.^2+2*z.^2 )*cos(theta0) );

dB_ix_dtheta = mu_0*M/4/pi*xyz_25.*( ( 2*x.^2-y.^2-z.^2 )*cos(theta0)*cos(phi0) +...
    3*x.*y*cos(theta0)*sin(phi0) -3*x.*z*sin(theta0) );

dB_iy_dtheta = mu_0*M/4/pi*xyz_25.*( 3*x.*y*cos(theta0)*cos(phi0) +...
    (-x.^2+2*y.^2-z.^2)*cos(theta0)*sin(phi0) - 3*y.*z*sin(theta0) );

dB_iz_dtheta = mu_0*M/4/pi*xyz_25.*( 3*x.*z*cos(theta0)*cos(phi0) + ...
    3*y.*z*cos(theta0)*sin(phi0) + ( x.^2+y.^2-2*z.^2 )*sin(theta0) );

dB_ix_dphi = mu_0*M/4/pi*xyz_25.*( 3*x.*y*sin(theta0)*cos(phi0) - ( 2*x.^2-y.^2-z.^2)*sin(theta0)*sin(phi0) );

dB_iy_dphi = mu_0*M/4/pi*xyz_25.*( (-x.^2+2*y.^2-z.^2)*sin(theta0)*cos(phi0) -3*x.*y*sin(theta0)*sin(phi0)  );

dB_iz_dphi = mu_0*M/4/pi*xyz_25.*(3*y.*z*sin(theta0)*cos(phi0) - 3*x.*z*sin(theta0)*sin(phi0) );

dB_ix_dm = mu_0/4/pi*xyz_25.*((2*x.^2-y.^2-z.^2)*sin(theta0)*cos(phi0) +3*x.*y*sin(theta0)*sin(phi0) + 3*x.*z*cos(phi0));

dB_iy_dm = mu_0/4/pi*xyz_25.*(3*x.*y*sin(theta0)*cos(phi0) + (-x.^2+2*y.^2-z.^2)*sin(theta0)*sin(phi0) + 3*y.*z*cos(phi0) );

dB_iz_dm = mu_0/4/pi*xyz_25.*(3*x.*z*sin(theta0)*cos(phi0) + 3*y.*z*sin(theta0)*sin(phi0) + (-x.^2-y.^2+2*z.^2)*cos(phi0) );

dB_mi_dx = sum(meas_dir_s.*[dB_ix_dx;dB_iy_dx;dB_iz_dx],1);
% dB_mi_dx = diag(meas_dir_s'*[dB_ix_dx;dB_iy_dx;dB_iz_dx]);

dB_mi_dy = sum(meas_dir_s.*[dB_ix_dy;dB_iy_dy;dB_iz_dy],1);
% dB_mi_dy = diag(meas_dir_s'*[dB_ix_dy;dB_iy_dy;dB_iz_dy]);

dB_mi_dz = sum(meas_dir_s.*[dB_ix_dz;dB_iy_dz;dB_iz_dz],1);
% dB_mi_dz = diag(meas_dir_s'*[dB_ix_dz;dB_iy_dz;dB_iz_dz]);

dB_mi_dtheta = sum(meas_dir_s.*[dB_ix_dtheta;dB_iy_dtheta;dB_iz_dtheta],1);
% dB_mi_dtheta = diag(meas_dir_s'*[dB_ix_dtheta;dB_iy_dtheta;dB_iz_dtheta]);

dB_mi_dphi = sum(meas_dir_s.*[dB_ix_dphi;dB_iy_dphi;dB_iz_dphi],1);
% dB_mi_dphi = diag(meas_dir_s'*[dB_ix_dphi;dB_iy_dphi;dB_iz_dphi]);

dB_mi_dm = sum(meas_dir_s.*[dB_ix_dm;dB_iy_dm;dB_iz_dm],1);
% dB_mi_dm = diag(meas_dir_s'*[dB_ix_dm;dB_iy_dm;dB_iz_dm]);

J = [ dB_mi_dx'  dB_mi_dy'  dB_mi_dz'  dB_mi_dtheta'  dB_mi_dphi' dB_mi_dm'];

end