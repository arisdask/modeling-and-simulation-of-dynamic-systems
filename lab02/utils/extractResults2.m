function [x, xdot, phi, theta, m_est, b_est, k_est, x_hat, xdot_hat, e_x, e_xdot] = extractResults2(t, y, lambda1, lambda2)
    % Extract states from ODE solver output for second-order filter version
    x = y(:, 1);                % Position
    xdot = y(:, 2);             % Velocity
    
    % Extract the actual phi values (not their derivatives)
    phi1 = y(:, 3);             % First filtered signal (without derivative)
    phi2 = y(:, 5);             % Second filtered signal (without derivative)
    phi3 = y(:, 7);             % Third filtered signal (without derivative)
    
    % Combine phi values for computation
    phi = [phi1, phi2, phi3];
    
    % Parameter estimates
    theta = y(:, 9:11);         % Parameter estimates - θ
    
    % Compute physical parameter estimates
    % Based on identification equation: y_hat = θ^T * φ = θ1*φ1 + θ2*φ2 + θ3*φ3
    % where θ1 = a1 - λ1 = b/m - λ1, θ2 = a2 - λ2 = k/m - λ2, θ3 = b0 = 1/m
    
    % From θ3 = 1/m, we get m = 1/θ3
    m_est = 1 ./ theta(:, 3);
    
    % From θ1 = b/m - λ1, we get b = m*(θ1 + λ1) = (θ1 + λ1)/θ3
    b_est = (theta(:, 1) + lambda1) ./ theta(:, 3);
    
    % From θ2 = k/m - λ2, we get k = m*(θ2 + λ2) = (θ2 + λ2)/θ3
    k_est = (theta(:, 2) + lambda2) ./ theta(:, 3);
    
    % Calculate predictions from model
    x_hat = sum(theta .* phi, 2);    % y_hat = θ^T * φ directly gives x_hat
    
    % Calculate xdot_hat by numerical differentiation of x_hat
    xdot_hat = gradient(x_hat, t);   % Compute velocity estimate by differentiating position estimate
    
    % Calculate estimation errors
    e_x = x - x_hat;
    e_xdot = xdot - xdot_hat;
end
