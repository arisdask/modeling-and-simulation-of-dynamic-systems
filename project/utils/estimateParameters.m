function [theta_est, x_hat] = estimateParameters(t, x_data, u_data, basis_func, lambda)
    % Real-time parameter estimation using gradient method
    
    N = length(t);
    
    % Ensure we have consistent time step
    if N > 1
        dt = t(2) - t(1);
    else
        dt = 0.01;  % Default time step
    end
    
    % Ensure data is in column vector format
    if size(x_data, 2) > size(x_data, 1)
        x_data = x_data';
    end
    if size(u_data, 2) > size(u_data, 1)
        u_data = u_data';
    end
    
    % Determine number of parameters from basis function
    phi_test   = basis_func(x_data(1), u_data(1));
    num_params = length(phi_test);
    
    theta = zeros(num_params, 1);
    
    % Learning rate
    gamma = 1;
    
    % Filter states
    x_filtered = 0;
    u_filtered = 0;
    
    x_hat = zeros(N, 1);
    
    % Performed without using ode
    for i = 1:N
        % Update filter states (first-order filter: dx_f/dt = -λ*x_f + input)
        if i == 1
            x_filtered = x_data(i);
            u_filtered = u_data(i);
        else
            x_filtered = x_filtered + dt * (-lambda * x_filtered + x_data(i));
            u_filtered = u_filtered + dt * (-lambda * u_filtered + u_data(i));
        end
        
        % Compute regressor vector
        phi = basis_func(x_filtered, u_filtered);
        
        % Ensure phi is column vector
        if size(phi, 2) > size(phi, 1)
            phi = phi';
        end
        
        y_filtered =  x_filtered;
        y_hat      =  phi' * theta;
        error      =  y_filtered - y_hat;
        
        % Update parameters using gradient descent
        theta = theta + gamma * error * phi * dt;

        x_hat(i) = y_hat;
    end
    
    theta_est = theta;
end
