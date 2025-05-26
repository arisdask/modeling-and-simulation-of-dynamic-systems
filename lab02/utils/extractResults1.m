function [x, xdot, phi, theta, m_est, b_est, k_est, x_hat, xdot_hat, e_x, e_xdot] = extractResults1(t, y, lambda1)
    % Extract states from ODE solver output
    x = y(:, 1);                % Position
    xdot = y(:, 2);             % Velocity
    phi = y(:, 3:5);            % Filtered signals
    theta = y(:, 6:8);          % Parameter estimates
    
    % Compute physical parameter estimates
    % Based on identification equation: y_hat = θ^T * φ = θ1*φ1 + θ2*φ2 + θ3*φ3
    % where θ1 = a1 - λ1 = b/m - λ1, θ2 = a2 = k/m, θ3 = b0 = 1/m
    
    % From θ3 = 1/m, we get m = 1/θ3
    m_est = 1 ./ theta(:, 3);
    
    % From θ1 = b/m - λ1, we get b = m*(θ1 + λ1) = (θ1 + λ1)/θ3
    b_est = (theta(:, 1) + lambda1) ./ theta(:, 3);
    
    % From θ2 = k/m, we get k = m*θ2 = θ2/θ3
    k_est = theta(:, 2) ./ theta(:, 3);
    
    % Calculate predictions from model
    xdot_hat = sum(theta .* phi, 2);    % y_hat = θ^T * φ
    
    % Integrate xdot_hat to get x_hat using MATLAB's cumtrapz function
    x_hat = cumtrapz(t, xdot_hat);      % Integrate velocity to get position
    
    % Alternative implementation without built-in function:
    % x_hat_alt = zeros(size(t));
    % for i = 2:length(t)
    %     % Simple trapezoidal integration
    %     dt = t(i) - t(i-1);
    %     x_hat_alt(i) = x_hat_alt(i-1) + 0.5 * (xdot_hat(i) + xdot_hat(i-1)) * dt;
    % end
    % x_hat = x_hat_alt;
    
    % Calculate estimation errors
    e_x = x - x_hat;                    % Not directly estimated
    e_xdot = xdot - xdot_hat;
end
