function dtheta_projected = applyProjection(theta, dtheta_unconstrained, Gamma, nabla_g1, nabla_g2, nabla_g3)
    % Apply projection method to handle constraints
    % Constraints: g1 = -a11 - 3 ≤ 0, g2 = a11 + 1 ≤ 0, g3 = -b2 + 1 ≤ 0
    
    % Evaluate constraint functions
    g1 = -theta(1) - 3;  % -a11 - 3
    g2 = theta(1) + 1;   % a11 + 1
    g3 = -theta(6) + 1;  % -b2 + 1
    
    % Tolerance for constraint boundary
    tol = 1e-2;
    
    % Determine active constraints
    active_constraints = [];
    
    % Check constraint g1 and g2
    if g1 > 0 || (abs(g1) <= tol && (dtheta_unconstrained' * nabla_g1) > 0)
        active_constraints = [active_constraints, 1];
    elseif g2 > 0 || (abs(g2) <= tol && (dtheta_unconstrained' * nabla_g2) > 0)
        active_constraints = [active_constraints, 2];
    end
    
    % Check constraint g3: b2 ≥ 1
    if g3 > 0 || (abs(g3) <= tol && (dtheta_unconstrained' * nabla_g3) > 0)
        active_constraints = [active_constraints, 3];
    end
    
    % Start with unconstrained update
    dtheta_projected = Gamma * dtheta_unconstrained;
    
    % Apply projection for each active constraint
    for i = active_constraints
        switch i
            case 1
                nabla_g = nabla_g1;
            case 2
                nabla_g = nabla_g2;
            case 3
                nabla_g = nabla_g3;
        end
        
        % Project out the component along the constraint gradient
        % Projection formula: dtheta_new = dtheta - Gamma * (nabla_g * nabla_g') / (nabla_g' * Gamma * nabla_g) * Gamma * dtheta
        denominator = nabla_g' * Gamma * nabla_g;
        
        if abs(denominator) > 1e-10  % Avoid division by zero
            projection_term = (Gamma * (nabla_g * nabla_g') * Gamma * dtheta_unconstrained) / denominator;
            dtheta_projected = dtheta_projected - projection_term;
        end
    end
end
