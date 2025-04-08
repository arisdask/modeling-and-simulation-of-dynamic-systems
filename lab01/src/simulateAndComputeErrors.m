function [t_est, x_est, e_q, e_q_dot] = simulateAndComputeErrors(pendulumODE, time_continuous, x0, m_est, L_est, c_est, g, inputTorque, t_true, x_true)
    % Simulates the system with estimated parameters and computes errors
    %
    % Inputs:
    %   pendulumODE - Function handle for the pendulum ODE
    %   time_continuous - Continuous time vector
    %   x0 - Initial conditions
    %   m_est, L_est, c_est - Estimated parameters
    %   g - Gravity constant
    %   inputTorque - Function handle for input torque
    %   t_true - True time vector
    %   x_true - True state trajectories
    %
    % Outputs:
    %   t_est - Time vector for the estimated trajectory
    %   x_est - Estimated state trajectories
    %   e_q - Error in angle
    %   e_q_dot - Error in angular velocity
    
    % Simulate the system with estimated parameters
    [t_est, x_est] = solveODE(pendulumODE, time_continuous, x0, m_est, L_est, c_est, g, inputTorque);
    
    % Interpolate to common timebase for comparison
    x_true_resampled = interp1(t_true, x_true, t_est, 'linear');
    
    % Compute errors
    e_q = x_true_resampled(:,1) - x_est(:,1);  % Error in q
    e_q_dot = x_true_resampled(:,2) - x_est(:,2);  % Error in q'
end
