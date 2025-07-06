function omega = computeBiasError(t, omega_bar, omega_freq)
    % Compute bounded bias error ω(t)
    % Model as combination of sinusoidal components with different frequencies
    
    omega1 = 0.8 * omega_bar * sin(omega_freq(1) * t);
    omega2 = 0.6 * omega_bar * cos(omega_freq(2) * t + pi/4);
    
    omega = [omega1; omega2];
    
    % Ensure ||ω(t)|| ≤ ω_bar
    omega_norm = norm(omega);
    if omega_norm > omega_bar
        omega = omega * (omega_bar / omega_norm);
    end
end
