function [pm] = PM_backward_estimation(y_meas,pm_init,xyz_s,meas_dir_s)

n_par = length(pm_init);

iteration = 0;  % current number of 
stop = 0;  % stop flag
MaxIter = 1000;   % maximum number of iteration
epsilon_1 = 1e-12;  % convergence tolerance for gradient
epsilon_2 = 1e-10;  % convergence tolerance for parameter
epsilon_3 = 1e-15;  % convergence tolerance for residual error
epsilon_4 = 1e-8;  % determines acceptance of a L-M step
pm = pm_init; % initialization on parameters
lambda_0 = 1e-3;  % user-specified initial lambda

% Jacobian initialization
J = Analytic_Jacobian_Mag_cir_array_test(pm_init,xyz_s,meas_dir_s);

% initialize L-M parameters
lambda = lambda_0*max(diag(J'*J));
v = 2;  % scale factor

% Main iteration starts 
% pm is the parameter evolving in each iteration
while ( ~stop && iteration <= MaxIter)
    iteration = iteration + 1;
    
    % evaluate function value using current parameter
    B_meas = PM_forward_field(pm,xyz_s,meas_dir_s);  % compute measured component
    B_meas = B_meas';  % reshape into column vector
    
    % compute residual error between model and measurement
    err = y_meas - B_meas;
    
    % check stopping criterion 3
    if ( err'*err < epsilon_3 &&  iteration > 2 )
        stop = 1;
        break
    end
    
    % compute gradient
    grad = J'*err;
    % check stopping criterion 1
    if ( max(abs(grad)) < epsilon_1  &&  iteration > 2 )
        stop = 1;
        break
    end
    
    % parameter increment change
    delta = (J'*J + lambda*eye(n_par))\J'*err;  % simplest version
    
    % check stopping criterion 2 
    if ( max(abs(delta)./(abs(pm)+1e-20)) < epsilon_2  &&  iteration > 2 )
        stop = 1;
        break
    end
    
    % check if [pm+delta] is better than [pm]
    pm_try = pm + delta;
    
    % evaluate value using pm_try
    B_meas_try = PM_forward_field(G,pm_try,xyz_s,meas_dir_s);  % compute measurement based on given direction
    B_meas_try = B_meas_try';  % reshape into column vector
    err_try = y_meas - B_meas_try;
    
    if ~all(isfinite(err_try))  % floating point error; break
        stop = 1;
        break
    end
    
    % compute rho for determining acceptance
    rho = (err'*err - err_try'*err_try)/(delta'*(lambda*delta + J'*err) );
    
    if rho > epsilon_4   % if better
        pm = pm_try;  % accpet pm_try
        lambda = lambda*max( 1/3, 1-(2*rho-1)^3 );
        v = 2;
        % update Jacobian matrix
        J = Analytic_Jacobian_Mag_cir_array_test(pm,xyz_s,meas_dir_s);
    else      % if not better
        lambda = lambda*v;
        v = 2*v;
    end
    % check stopping criterion 4
    if ( iteration == MaxIter )
        disp(' !! Maximum Number of Iterations Reached Without Convergence !!')
        stop = 1;
    end
    

end
end