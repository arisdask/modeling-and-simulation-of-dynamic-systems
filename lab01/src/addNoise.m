function [q_noisy, q_dot_noisy, u_sampled_noisy] = addNoise(q, q_dot, u_sampled, noise_level_q, noise_level_q_dot)
    % Adds white Gaussian noise to the angle and angular velocity measurements
    %
    % Inputs:
    %   q - Clean angle measurements
    %   q_dot - Clean angular velocity measurements
    %   u_sampled - Clean input measurements
    %   noise_level_q - Standard deviation of noise for angle
    %   noise_level_q_dot - Standard deviation of noise for angular velocity
    
    % Generate white Gaussian noise
    noise_q = noise_level_q * randn(size(q));
    noise_q_dot = noise_level_q_dot * randn(size(q_dot));
    noise_u = noise_level_q * randn(size(u_sampled));
    
    % Add noise to measurements
    q_noisy = q + noise_q;
    q_dot_noisy = q_dot + noise_q_dot;
    u_sampled_noisy = u_sampled + noise_u;
end
