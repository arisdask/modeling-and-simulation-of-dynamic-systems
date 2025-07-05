function dxdt = biasedSystemDynamics(t, x, A_true, B_true, lambda, Gamma, M, sigma, omega_bar, omega_freq)
    % Implements the gradient method with σ-modification for biased systems
    %
    % State vector x = [x_true; theta_lambda; filter_states]
    % where:
    % - x_true: true system states (2×1)
    % - theta_lambda: filtered parameter estimates (6×1)
    % - filter_states: internal states of all filters (3×1)
    
    % Extract states
    x_true = x(1:2);
    theta_lambda = x(3:8);
    filter_states = x(9:11);
    
    % Input signal u(t) - rich enough signal for parameter estimation
    u = 2.5 * sin(t) + 1.5 * cos(1.5*t) + 0.8 * sin(3*t);
    
    % Bias error ω(t) - bounded disturbance
    omega = computeBiasError(t, omega_bar, omega_freq);
    
    % True system dynamics with bias: dx/dt = A_true*x + B_true*u + ω
    dx_true_dt = A_true * x_true + B_true * u + omega;
    
    % Extract filtered signals from filter states
    x1_filt = filter_states(1);
    x2_filt = filter_states(2);
    u_filt  = filter_states(3);
    
    % Construct regressor matrix Φ(t) using filtered signals
    Phi = [x1_filt, x2_filt, 0,       0,       u_filt, 0;
           0,       0,       x1_filt, x2_filt, 0,      u_filt];
    
    % Compute system output y = [x1; x2]
    y = x_true;
    
    % Compute estimated output: ŷ = Φ * θ_λ
    y_hat = Phi * theta_lambda;
    
    % Compute estimation error: e = y - ŷ
    e = y - y_hat;
    
    % Compute gradient of cost function: ∇K = -Φ^T * e
    grad_K = -Phi' * e;
    
    % Compute σ-modification term
    theta_norm = norm(theta_lambda);
    sigma_delta = computeSigmaModification(theta_norm, M, sigma);
    
    % Parameter update law with σ-modification:
    % dθ_λ/dt = -σ_δ * Γ * θ_λ - Γ * ∇K
    dtheta_lambda_dt = -sigma_delta * Gamma * theta_lambda - Gamma * grad_K;
    
    % Apply constraints through projection
    dtheta_lambda_dt = applyProjectionThetaLambda(theta_lambda, dtheta_lambda_dt, lambda);
    
    % Update filter states
    dfilter_states_dt = updateFilterStates(filter_states, x_true, u, lambda);
    
    % Combine all derivatives
    dxdt = [dx_true_dt; dtheta_lambda_dt; dfilter_states_dt];
end
