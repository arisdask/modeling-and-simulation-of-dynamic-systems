function printRMSErrors(e_q, e_q_dot)
    % Calculates and prints RMS errors
    %
    % Inputs:
    %   e_q - Error in angle
    %   e_q_dot - Error in angular velocity
    
    rms_e_q = sqrt(mean(e_q.^2));
    rms_e_q_dot = sqrt(mean(e_q_dot.^2));
    
    fprintf('\nRMS Error:\n');
    fprintf('RMS e_q = %.6f, RMS e_q_dot = %.6f\n', rms_e_q, rms_e_q_dot);
end
