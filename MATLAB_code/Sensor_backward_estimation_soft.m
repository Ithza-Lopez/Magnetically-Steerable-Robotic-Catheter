function [xyz_s_opt,O_opt,M_opt] = Sensor_backward_estimation_soft(B_meas,pm_all,xyz_s_0,O_0,meas_dir_s,M_0)
% use Levenberg–Marquardt algorithm: minimize squared error
% weighting among signals needed since x,y,z are in [m] unit but theta,phi
% are in [rad] unit

% regularization TBD

% input:
% xyz_s_0: sensor configuration guess <- 3x1
% O: 3x3 matrix initial guess
% pm_all: all permanent magnet positions&orientations <- 5xN 
% y_meas: measurement data of one sensor for all PM positions <- (3xN)x1 vector 
% M_0: magnetic moment
% meas_dir_s: sensor measurement direction information <- R^(3 x 3)

% output:
% xyz_s_opt : optimal sensor position
% O_opt: optimal O matrix
% M_opt: optimal magnetic moment
scale = 1;
O_0_vec = reshape(O_0,[1,9])'/scale;

n_par = length(xyz_s_0)+length(O_0_vec);

% form the parameter vector: [x_s;y_s;z_s;O11;O12;O13;O21;O22;O23;O31;O32;O33]
para_0 = [xyz_s_0;O_0_vec];

iteration = 0;  % current number of iterations
stop = 0;  % stop flag
MaxIter = 500;   % maximum number of iteration
epsilon_1 = 1e-12;  % convergence tolerance for gradient
epsilon_2 = 1e-12;  % convergence tolerance for parameter
epsilon_3 = 1e-16;  % convergence tolerance for residual error
epsilon_4 = 1e-8;  % determines acceptance of a L-M step
para = para_0; % initialization on parameters
lambda_0 = 1e-3;  % user-specified initial lambda
M = M_0;
xyz_s = xyz_s_0; 
O = O_0/scale;

% Jacobian initialization
J = Analytic_Jacobian_Sensor_cir_array_soft(pm_all,xyz_s_0,meas_dir_s,M_0,B_meas);

% initialize L-M parameters
lambda = lambda_0*max(diag(J'*J));
v = 2;  % scale factor

% Main iteration starts 
% xyz is the parameter evolving in each iteration
while ( ~stop && iteration <= MaxIter)
    iteration = iteration + 1;
    
    % evaluate function value using current parameter
    B_meas_sim = Sensor_forward(pm_all,xyz_s,meas_dir_s,M);  % compute measured component
    
    % compute residual error between model and measurement
    err = mod_B_meas_sensor(O*scale,B_meas) - B_meas_sim;
    
    % check stopping criterion 3
    if ( err'*err < epsilon_3 &&  iteration > 2 )
        fprintf('Convergence in residual error')
        stop = 1;
        break
    end
    
    % compute gradient
    grad = J'*err;
    
    % check stopping criterion 1
    if ( max(abs(grad)) < epsilon_1  &&  iteration > 2 )
        fprintf(' Convergence in gradient')
        stop = 1;
        break
    end
    
    % parameter increment change
    delta = (J'*J + lambda*eye(n_par))\J'*err;  % simplest version
    
    % check stopping criterion 2 
    
    if ((max(abs(delta)./(abs(para)+1e-20))) < epsilon_2  &&  iteration > 2 )
        fprintf('Convergence in Parameters')
        stop = 1;
        break
    end
    
    % check if [pm+delta] is better than [pm]
    para_try = para + delta;
    xyz_s_try = para_try(1:3);
    O_vec_try = para_try(4:end);
    O_try = [O_vec_try(1:3)';O_vec_try(4:6)';O_vec_try(7:9)'];
    
    % evaluate value using pm_try
    B_meas_try = Sensor_forward(pm_all,xyz_s_try,meas_dir_s,M);  % compute measurement based on given direction
    err_try = mod_B_meas_sensor(O_try*scale,B_meas) - B_meas_try;
    
    if ~all(isfinite(err_try))  % floating point error; break
        stop = 1;
        break
    end
    
    % compute rho for determining acceptance
    rho = (err'*err - err_try'*err_try)/(delta'*(lambda*delta + J'*err) );
    
    if rho > epsilon_4   % if better
        para = para_try;  % accpet pm_try
        xyz_s = para_try(1:3);
        O = O_try;
        lambda = lambda*max( 1/3, 1-(2*rho-1)^3 );
        v = 2;
        % update Jacobian matrix
        J = Analytic_Jacobian_Sensor_cir_array_soft(pm_all,xyz_s,meas_dir_s,M,B_meas);
    else      % if not better
        lambda = lambda*v;
        v = 2*v;
    end
    
    % print history
    fprintf('%6s %9s %9s %9s\n',...
        'iter', '||error||^2', '||grad||', 'lambda');
    fprintf('%6i %9.2e %9.2e %9.2e\n',...
        iteration, err'*err, norm(grad), lambda);
    
    % check stopping criterion 4
    if ( iteration == MaxIter )
        disp(' !! Maximum Number of Iterations Reached Without Convergence !!')
        stop = 1;
    end
    

end
xyz_s_opt = para(1:3);
O_vec_opt = para(4:end)*scale;
O_opt = [O_vec_opt(1:3)';O_vec_opt(4:6)';O_vec_opt(7:9)'];
M_opt = M;
end