function [x, xdot, x_hat, xdot_hat, m_est, b_est, k_est, e_x, x_measured] = extractResultsWithNoise(y, t, n0, f0)
    % Extract states
    x = y(:, 1);           % True position
    xdot = y(:, 2);        % True velocity
    x_hat = y(:, 3);       % Estimated position
    xdot_hat = y(:, 4);    % Estimated velocity
    theta = y(:, 5:7);     % Parameter estimates
    
    % Add noise to position measurement
    noise = n0 * sin(2*pi*f0*t);
    x_measured = x + noise';
    
    % Compute physical parameters
    m_est = 1 ./ theta(:, 3);
    b_est = theta(:, 1) .* m_est;
    k_est = theta(:, 2) .* m_est;
    
    % Compute estimation error using noisy measurement
    e_x = x_measured - x_hat;
end