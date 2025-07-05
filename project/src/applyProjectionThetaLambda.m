function dtheta_projected = applyProjectionThetaLambda(theta_lambda, dtheta_unconstrained, lambda)
    % Apply projection to handle constraints on original parameters
    % Converts θ_λ back to original parameters for constraint checking !!
    
    % Extract original parameters from filtered parameters
    a11 = theta_lambda(1) - lambda;  % a11 = (a11 + λ) - λ
    b2 = theta_lambda(6);            % b2 unchanged
    
    % Constraint functions:
    % g1 = -a11 - 3 ≤ 0  (a11 ≥ -3)
    % g2 = a11 + 1 ≤ 0   (a11 ≤ -1)
    % g3 = -b2 + 1 ≤ 0   (b2 ≥ 1)
    
    g1 = -a11 - 3;
    g2 = a11 + 1;
    g3 = -b2 + 1;
    
    % Constraint gradients in terms of θ_λ (same as in 1a)
    nabla_g1 = [1; 0; 0; 0; 0; 0];
    nabla_g2 = [-1; 0; 0; 0; 0; 0];
    nabla_g3 = [0; 0; 0; 0; 0; -1];
    
    % Start with unconstrained update
    dtheta_projected = dtheta_unconstrained;
    
    % Tolerance for constraint boundary
    tol = 1e-3;
    
    % Apply projection for active constraints
    if g1 > 0 || (abs(g1) <= tol && (dtheta_unconstrained' * nabla_g1) > 0)
        projection = ((nabla_g1' * dtheta_unconstrained) / (nabla_g1' * nabla_g1)) * nabla_g1;
        dtheta_projected = dtheta_projected - projection;
    elseif g2 > 0 || (abs(g2) <= tol && (dtheta_unconstrained' * nabla_g2) > 0)
        projection = ((nabla_g2' * dtheta_unconstrained) / (nabla_g2' * nabla_g2)) * nabla_g2;
        dtheta_projected = dtheta_projected - projection;
    end
    
    if g3 > 0 || (abs(g3) <= tol && (dtheta_unconstrained' * nabla_g3) > 0)
        projection = ((nabla_g3' * dtheta_unconstrained) / (nabla_g3' * nabla_g3)) * nabla_g3;
        dtheta_projected = dtheta_projected - projection;
    end
end
