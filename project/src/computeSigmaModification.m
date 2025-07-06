function sigma_delta = computeSigmaModification(theta_norm, M, sigma)
    % Compute σ-modification parameter based on parameter norm
    % σ_δ = 0 if ||θ|| < M
    % σ_δ = σ(||θ||/M - 1) if M <= ||θ|| <= 2M  
    % σ_δ = σ if ||θ|| > 2M
    
    if theta_norm < M
        sigma_delta = 0;
    elseif theta_norm <= 2*M
        sigma_delta = sigma * (theta_norm/M - 1);
    else
        sigma_delta = sigma;
    end
end
