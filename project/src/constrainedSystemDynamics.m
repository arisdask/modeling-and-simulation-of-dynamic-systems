function dxdt = constrainedSystemDynamics(t, x, A_true, B_true, C, Gamma, nabla_g1, nabla_g2, nabla_g3)
    % Extract states from the input vector
    % x = [x1; x2; x_hat1; x_hat2; theta1; theta2; theta3; theta4; theta5; theta6]
    
    % True system states
    x_true = x(1:2);
    
    % Estimated states
    x_hat = x(3:4);
    
    % Parameter estimates: theta = [a11, a12, a21, a22, b1, b2]'
    theta = x(5:10);
    
    % Input signal u(t) - rich enough signal for parameter estimation
    u = 2.5 * sin(t) + 1.5 * cos(1.5*t) + 0.8 * sin(3*t);
    
    % Reconstruct estimated matrices from parameters
    A_hat = [theta(1), theta(2);
             theta(3), theta(4)];
    
    B_hat = [theta(5); theta(6)];
    
    % Compute estimation error
    e = x_true - x_hat;
    
    % True system dynamics: dx/dt = A_true * x + B_true * u
    dx_true_dt = A_true * x_true + B_true * u;
    
    % Mixed structure estimator dynamics: 
    % dx_hat/dt = A_hat * x_true + B_hat * u + C * e
    dx_hat_dt = A_hat * x_true + B_hat * u + C * e;
    
    % Compute unconstrained parameter update
    % Based on: dtheta/dt = [e1*x1; e1*x2; e2*x1; e2*x2; e1*u; e2*u]
    dtheta_unconstrained = [e(1) * x_true(1);  % da11/dt
                           e(1) * x_true(2);  % da12/dt  
                           e(2) * x_true(1);  % da21/dt
                           e(2) * x_true(2);  % da22/dt
                           e(1) * u;          % db1/dt
                           e(2) * u];         % db2/dt
    
    % Apply projection to handle constraints
    dtheta_dt = applyProjection(theta, dtheta_unconstrained, Gamma, nabla_g1, nabla_g2, nabla_g3);
    
    % Combine all derivatives
    dxdt = [dx_true_dt;   % True system dynamics
            dx_hat_dt;    % Estimator dynamics
            dtheta_dt];   % Parameter adaptation with constraints
end
