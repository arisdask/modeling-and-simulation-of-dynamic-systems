function printEstimationResults(param_name, est_value, true_value)
    % Prints the comparison between estimated and true parameter values
    % 
    % Inputs:
    %   param_name - Name of the parameter (string)
    %   est_value - Estimated parameter value
    %   true_value - True parameter value
    
    error_percent = 100 * abs(est_value - true_value) / true_value;
    fprintf('%s_est = %.4f (true: %.4f, error: %.2f%%)\n', param_name, est_value, true_value, error_percent);
end
