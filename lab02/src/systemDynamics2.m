%% Function to define system dynamics
function dydt = systemDynamics2(t, y, m_true, b_true, k_true, lambda1, lambda2, gamma, input_type)
    % Extract states from the input vector
    x = y(1);          % Position - x
    xdot = y(2);       % Velocity - x'
    phi1 = y(3:4);     % First filtered signal (phi1 and its derivative)
    phi2 = y(5:6);     % Second filtered signal (phi2 and its derivative)
    phi3 = y(7:8);     % Third filtered signal (phi3 and its derivative)
    theta = y(9:11);   % Parameter estimates - θ
    
    % Determine input based on the specified case
    if input_type == 1
        u = 2.5;          % Constant input
    else
        u = 2.5 * sin(t); % Sinusoidal input
    end
    
    % System dynamics: d²x/dt² = (u - b*dx/dt - k*x) / m
    acceleration = (u - b_true*xdot - k_true*x) / m_true;
    
    % Second-order filter dynamics for each regressor:
    % For φ₁
    phi1_dot = phi1(2);  % Extract the derivative from state vector
    phi1_ddot = -lambda1*phi1_dot - lambda2*phi1(1) - xdot;  % φ₁: φ̈₁ = -λ₁φ̇₁ - λ₂φ₁ - ẋ
    
    % For φ₂
    phi2_dot = phi2(2);  % Extract the derivative from state vector
    phi2_ddot = -lambda1*phi2_dot - lambda2*phi2(1) - x;  % φ₂: φ̈₂ = -λ₁φ̇₂ - λ₂φ₂ - x
    
    % For φ₃
    phi3_dot = phi3(2);  % Extract the derivative from state vector
    phi3_ddot = -lambda1*phi3_dot - lambda2*phi3(1) + u;  % φ₃: φ̈₃ = -λ₁φ̇₃ - λ₂φ₃ + u
    
    % Combine regressors into a single vector for output calculation
    phi = [phi1(1); phi2(1); phi3(1)];
    
    % Identifier output and error
    y_hat = theta' * phi;
    e = x - y_hat;  % e = y - yHat (where y is the x(t) in this analysis)
    
    % Gradient update for parameter estimation
    dtheta = gamma * e * phi;
    
    % Combine all derivatives into a single vector
    dydt = [xdot; 
            acceleration; 
            phi1_dot; phi1_ddot; 
            phi2_dot; phi2_ddot; 
            phi3_dot; phi3_ddot; 
            dtheta];
end
