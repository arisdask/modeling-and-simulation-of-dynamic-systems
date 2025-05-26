%% Function to define system dynamics
function dydt = systemDynamics1(t, y, m_true, b_true, k_true, lambda1, gamma, input_type)
    % Extract states from the input vector
    x = y(1);          % Position - x
    xdot = y(2);       % Velocity - x'
    phi = y(3:5);      % Filtered signals - φ
    theta = y(6:8);    % Parameter estimates - θ
    
    % Determine input based on the specified case
    if input_type == 1
        u = 2.5;          % Constant input
    else
        u = 2.5 * sin(t); % Sinusoidal input
    end
    
    % System dynamics: d²x/dt² = (u - b*dx/dt - k*x) / m
    acceleration = (u - b_true*xdot - k_true*x) / m_true;
    
    % Filter dynamics: dphi/dt = -lambda1*phi + [-xdot; -x; u]
    dphi = zeros(3, 1);
    dphi(1) = -lambda1*phi(1) - xdot;
    dphi(2) = -lambda1*phi(2) - x;
    dphi(3) = -lambda1*phi(3) + u;
    
    % Identifier output and error
    y_hat = theta' * phi;
    e = xdot - y_hat;  % e = y - yHat (where y is the x'(t) in our analysis)
    
    % Gradient update for parameter estimation
    dtheta = gamma * e * phi;
    
    % Combine all derivatives into a single vector
    dydt = [xdot; acceleration; dphi; dtheta];
end
