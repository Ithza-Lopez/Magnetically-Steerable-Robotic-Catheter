function [xyz_s,M] = Sensor_backward_estimation_noM(y_meas,pm_all,xyz_s_0,meas_dir_s,M_0)
% use Levenberg–Marquardt algorithm: minimize squared error
% weighting among signals needed since x,y,z are in [m] unit but theta,phi
% are in [rad] unit

% regularization TBD

% input:
% xyz_s_0: sensor configuration guess <- 3x1
% pm_all: all permanent magnet positions&orientations <- 5xN 
% y_meas: measurement data of one sensor for all PM positions <- (3xN)x1 vector 
% M_0: magnetic moment
% meas_dir_s: sensor measurement direction information <- R^(3 x 3)

% output:
% xyz_s : optimal sensor position
% M: optimal magnetic moment

n_par = length(xyz_s_0);

iteration = 0;  % current number of iterations
stop = 0;  % stop flag
MaxIter = 100;   % maximum number of iteration
epsilon_1 = 1e-18;  % convergence tolerance for gradient
epsilon_2 = 1e-20;  % convergence tolerance for parameter
epsilon_3 = 1e-17;  % convergence tolerance for residual error
epsilon_4 = 1e-8;  % determines acceptance of a L-M step
xyz_s = xyz_s_0; % initialization on parameters
M = M_0; % initialization on parameters
lambda_0 = 1e-3;  % user-specified initial lambda

% Jacobian initialization
J = Analytic_Jacobian_Sensor_cir_array_noM(pm_all,xyz_s_0,meas_dir_s,M_0);

% initialize L-M parameters
lambda = lambda_0*max(diag(J'*J));
v = 2;  % scale factor

% Main iteration starts 
% xyz is the parameter evolving in each iteration
while ( ~stop && iteration <= MaxIter)
    iteration = iteration + 1;
    
    % evaluate function value using current parameter
    B_meas = Sensor_forward(pm_all,xyz_s,meas_dir_s,M);  % compute measured component
    
    % compute residual error between model and measurement
    err = y_meas - B_meas;
    
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
    
    if ((max(abs(delta)./(abs(xyz_s)+1e-20))) < epsilon_2  &&  iteration > 2 )
        fprintf('Convergence in Parameters')
        stop = 1;
        break
    end
    
    % check if [pm+delta] is better than [pm]
    xyz_s_try = xyz_s + delta;
    
    % evaluate value using pm_try
    B_meas_try = Sensor_forward(pm_all,xyz_s_try,meas_dir_s,M);  % compute measurement based on given direction
    err_try = y_meas - B_meas_try;
    
    if ~all(isfinite(err_try))  % floating point error; break
        stop = 1;
        break
    end
    
    % compute rho for determining acceptance
    rho = (err'*err - err_try'*err_try)/(delta'*(lambda*delta + J'*err) );
    
    if rho > epsilon_4   % if better
        xyz_s = xyz_s_try;  % accpet pm_try
        lambda = lambda*max( 1/3, 1-(2*rho-1)^3 );
        v = 2;
        % update Jacobian matrix
        J = Analytic_Jacobian_Sensor_cir_array_noM(pm_all,xyz_s,meas_dir_s,M);
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
end