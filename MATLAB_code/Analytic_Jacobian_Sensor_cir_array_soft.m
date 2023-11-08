function [J] = Analytic_Jacobian_Sensor_cir_array_soft(pm0,xyz_s,meas_dir_s,M,B_meas)

% input
% pm0 : 5x(3xN) matrix magnet setting at each slot location
% xyz_s: 3x1 initial guess sensor location
% mear_dir_s: normalized sensor measurement direction information <- R^(3 x 3)
% M: current magnetic moment, assume using only one type of magnet
% B_meas: measurement vector <- R^(3N x 1 )

% output:
%        J : Jacobian w.r.t Sx,Sy,Sz  <- R^(3N x (3+9))

mu_0 = 4*pi*1e-7;
[~,N] = size(pm0);  % number of permanent magnet positions
[~,N_d] = size(meas_dir_s);  % count number of direction
p =  xyz_s - pm0(1:3,:) ; % position vector for all sensors
theta0 = pm0(4,:);  % current theta
phi0 = pm0(5,:);  % current phi

% position components at i-th sensor x_ij bar
x = p(1,:);
y = p(2,:);
z = p(3,:);

% pre-computations
xyz_sqr = sum(p.^2,1); %(x_i^2 + y_i^2 + z_i^2)
xyz_35 = xyz_sqr.^-3.5;
xyz_25 = xyz_sqr.^-2.5;


% matrix creation
dB_ix_dx = -mu_0*M/4/pi*xyz_35.*( 3*x.*( 2*x.^2-3*y.^2-3*z.^2 ).*sin(theta0).*cos(phi0) + ...
    3*y.*( 4*x.^2-y.^2-z.^2 ).*sin(theta0).*sin(phi0) + 3*z.*( 4*x.^2-y.^2-z.^2 ).*cos(theta0) );

dB_iy_dx = -mu_0*M/4/pi*xyz_35.*( 3*y.*( 4*x.^2-y.^2-z.^2 ).*sin(theta0).*cos(phi0) + ...
    3*x.*( -x.^2+4*y.^2-z.^2 ).*sin(theta0).*sin(phi0) + 15*x.*y.*z.*cos(theta0) );

dB_iz_dx = -mu_0*M/4/pi*xyz_35.*( 3*z.*( 4*x.^2-y.^2-z.^2 ).*sin(theta0).*cos(phi0) + ...
    15*x.*y.*z.*sin(theta0).*sin(phi0) + 3*x.*( -x.^2-y.^2+4*z.^2 ).*cos(theta0) );

dB_ix_dy = -mu_0*M/4/pi*xyz_35.*( 3*y.*( 4*x.^2-y.^2-z.^2 ).*sin(theta0).*cos(phi0) + ...
    3*x.*( -x.^2+4*y.^2-z.^2 ).*sin(theta0).*sin(phi0) + 15*x.*y.*z.*cos(theta0) );

dB_iy_dy = -mu_0*M/4/pi*xyz_35.*( 3*x.*( -x.^2+4*y.^2-z.^2 ).*sin(theta0).*cos(phi0) + ...
    3*y.*( -3*x.^2+2*y.^2-3*z.^2 ).*sin(theta0).*sin(phi0) + 3*z.*( -x.^2+4*y.^2-z.^2 ).*cos(theta0)   );

dB_iz_dy = -mu_0*M/4/pi*xyz_35.*( 15*x.*y.*z.*sin(theta0).*cos(phi0) + ...
    3*z.*( -x.^2+4*y.^2-z.^2 ).*sin(theta0).*sin(phi0) + 3*y.*( -x.^2-y.^2+4*z.^2 ).*cos(theta0) );

dB_ix_dz = -mu_0*M/4/pi*xyz_35.*( 3*z.*( 4*x.^2-y.^2-z.^2 ).*sin(theta0).*cos(phi0) + ...
    15*x.*y.*z.*sin(theta0).*sin(phi0) + 3*x.*( -x.^2-y.^2+4*z.^2 ).*cos(theta0) );

dB_iy_dz = -mu_0*M/4/pi*xyz_35.*( 15*x.*y.*z.*sin(theta0).*cos(phi0) + ...
    3*z.*( -x.^2+4*y.^2-z.^2 ).*sin(theta0).*sin(phi0) + 3*y.*( -x.^2-y.^2+4*z.^2 ).*cos(theta0)  );

dB_iz_dz = -mu_0*M/4/pi*xyz_35.*( 3*x.*(-x.^2-y.^2+4*z.^2).*sin(theta0).*cos(phi0) + ...
    3*y.*( -x.^2-y.^2+4*z.^2 ).*sin(theta0).*sin(phi0) + 3*z.*( -3*x.^2-3*y.^2+2*z.^2 ).*cos(theta0) );

dB_mi_dx = reshape(meas_dir_s'*[dB_ix_dx;dB_iy_dx;dB_iz_dx],[N_d*N,1]);

dB_mi_dy = reshape(meas_dir_s'*[dB_ix_dy;dB_iy_dy;dB_iz_dy],[N_d*N,1]);

dB_mi_dz = reshape(meas_dir_s'*[dB_ix_dz;dB_iy_dz;dB_iz_dz],[N_d*N,1]);

temp = reshape(B_meas,[3,N])';

dB_meas_Os = zeros(3*N,9);

for i = 1:N
    dB_meas_Os((i*3-2),1:3) = -temp(i,:);
    dB_meas_Os((i*3-1),4:6) = -temp(i,:);
    dB_meas_Os((i*3),7:9) = -temp(i,:);
end 


J = [ dB_mi_dx  dB_mi_dy  dB_mi_dz dB_meas_Os];

end