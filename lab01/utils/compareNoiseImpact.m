function compareNoiseImpact(param_name, clean_value, noisy_value, true_value)
    % Compares parameter estimates with and without noise
    %
    % Inputs:
    %   param_name - Name of the parameter (string)
    %   clean_value - Parameter estimate without noise
    %   noisy_value - Parameter estimate with noise
    %   true_value - True parameter value
    
    error_clean = 100 * abs(clean_value - true_value) / true_value;
    error_noisy = 100 * abs(noisy_value - true_value) / true_value;
    
    fprintf('%s: true=%.4f, clean=%.4f (error: %.2f%%), noisy=%.4f (error: %.2f%%)\n', ...
        param_name, true_value, clean_value, error_clean, noisy_value, error_noisy);
end
