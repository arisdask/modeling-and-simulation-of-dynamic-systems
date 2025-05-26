function dydt = parallelStructureDynamics(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3)
    % Extract states from the input vector
    x = y(1);          % True position
    xdot = y(2);       % True velocity
    x_hat = y(3);      % Estimated position
    xdot_hat = y(4);   % Estimated velocity
    theta = y(5:7);    % Parameter estimates [θ1 θ2 θ3]
    
    % Input u(t) = 2.5sin(t)
    u = 2.5 * sin(t);
    
    % True system dynamics: d²x/dt² = (u - b*dx/dt - k*x) / m
    acceleration = (u - b_true*xdot - k_true*x) / m_true;
    
    % Compute the estimation error and its derivative
    % e = x - x_hat;
    edot = xdot - xdot_hat;
    
    % Estimated system dynamics: d²x̂/dt² = -θ1*dx̂/dt - θ2*x̂ + θ3*u
    acceleration_hat = -theta(1)*xdot_hat - theta(2)*x_hat + theta(3)*u;
    
    % Parameter adaptation laws
    % dθ1/dt = -γ1*ė*dx̂/dt
    % dθ2/dt = -γ2*ė*x̂
    % dθ3/dt = +γ3*ė*u
    dtheta = zeros(3,1);
    dtheta(1) = -gamma1 * edot * xdot_hat;
    dtheta(2) = -gamma2 * edot * x_hat;
    dtheta(3) = gamma3 * edot * u;
    
    % Combine all derivatives into a single vector
    dydt = [xdot;                  % dx/dt
            acceleration;          % d²x/dt²
            xdot_hat;             % dx̂/dt
            acceleration_hat;      % d²x̂/dt²
            dtheta];              % dθ/dt
end
