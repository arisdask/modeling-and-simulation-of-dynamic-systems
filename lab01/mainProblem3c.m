% Modeling and Simulation of Dynamic Systems - Assignment 1
% Estimation of Unknown Parameters - Least Squares Method
% Author: Aristeidis Daskalopoulos - AEM:10640

%% Problem 3c: Effect of varying input amplitude on parameter estimation accuracy
% Initialize Problem
clear; clc; close all;

addpath('src');
addpath('utils');

% True system parameters
m_true = 0.75; L_true = 1.25; c_true = 0.15;
g = 9.81; omega = 2;

% Fixed sampling period
Ts = 0.1; % [sec]

% Range of input amplitudes to test
input_amplitudes = 1:1:10;  % Varying A0 from 1 to 10
num_amplitudes = length(input_amplitudes);

T_final = 20; % total simulation time [sec]
dt = 1e-4; % integration step

% Initial conditions
x0 = [0; 0]; % [q(0) q'(0)] = [0 0]

% Initialize arrays to store parameter estimation errors
m_error_2a = zeros(num_amplitudes, 1);
L_error_2a = zeros(num_amplitudes, 1);
c_error_2a = zeros(num_amplitudes, 1);
m_error_2b = zeros(num_amplitudes, 1);
L_error_2b = zeros(num_amplitudes, 1);
c_error_2b = zeros(num_amplitudes, 1);

% Define filter parameters
lambda1_2a = 0.2;     % For method 2a
lambda1_2b = 0.6;     % For method 2b
lambda2_2b = 0.08;    % For method 2b

% Loop through each input amplitude
for i = 1:num_amplitudes
    A0 = input_amplitudes(i);
    
    % Input torque function with current amplitude
    inputTorque = @(t) A0 * sin(omega * t);
    
    % Pendulum ODE function
    pendulumODE = @(t, x, m, L, c, g, inputFunc) [
        x(2); % q' = dq/dt = x_1'
        (inputFunc(t) - c*x(2) - m*g*L*x(1))/(m*L^2) % q" = x_2'
    ];
    
    % Generate "true" system data for current amplitude
    [t_true, x_true] = solveODE(pendulumODE, 0:dt:T_final, x0, m_true, L_true, c_true, g, inputTorque);
    
    % Sample the data at fixed Ts intervals
    sampled_indices = 1:round(Ts/dt):length(t_true);
    t_sampled = t_true(sampled_indices);
    x_sampled = x_true(sampled_indices, :);
    q_sampled = x_sampled(:, 1);
    q_dot_sampled = x_sampled(:, 2);
    u_sampled = arrayfun(inputTorque, t_sampled);
    
    %% Method 2a: Both q(t) and q'(t) are measurable
    % Create filter transfer functions
    s = tf('s');
    Lambda_s_2a = s + lambda1_2a;
    
    % Create filtered signals
    zeta1_2a = lsim(-1/Lambda_s_2a, q_dot_sampled, t_sampled);
    zeta2_2a = lsim(-1/Lambda_s_2a, q_sampled, t_sampled);
    zeta3_2a = lsim(1/Lambda_s_2a, u_sampled, t_sampled);
    
    % Create regression matrix and apply least squares
    zeta_matrix_2a = [zeta1_2a, zeta2_2a, zeta3_2a];
    theta_lambda_2a = (zeta_matrix_2a' * zeta_matrix_2a) \ (zeta_matrix_2a' * q_dot_sampled);
    
    % Extract parameters
    a1_minus_lambda1_2a = theta_lambda_2a(1);
    a2_2a = theta_lambda_2a(2);
    b0_2a = theta_lambda_2a(3);
    a1_2a = a1_minus_lambda1_2a + lambda1_2a;
    
    % Compute physical parameters
    c_est_2a = a1_2a/b0_2a;
    L_est_2a = g / a2_2a;
    m_est_2a = 1 / (b0_2a * L_est_2a^2);
    
    % Compute relative errors
    m_error_2a(i) = abs(m_est_2a - m_true) / m_true * 100;
    L_error_2a(i) = abs(L_est_2a - L_true) / L_true * 100;
    c_error_2a(i) = abs(c_est_2a - c_true) / c_true * 100;
    
    %% Method 2b: Only q(t) is measurable
    % Create filter transfer functions
    Lambda_s_2b = s^2 + lambda1_2b*s + lambda2_2b;
    
    % Create filtered signals
    zeta1_2b = lsim(-s/Lambda_s_2b, q_sampled, t_sampled);
    zeta2_2b = lsim(-1/Lambda_s_2b, q_sampled, t_sampled);
    zeta3_2b = lsim(1/Lambda_s_2b, u_sampled, t_sampled);
    
    % Create regression matrix and apply least squares
    zeta_matrix_2b = [zeta1_2b, zeta2_2b, zeta3_2b];
    theta_lambda_2b = (zeta_matrix_2b' * zeta_matrix_2b) \ (zeta_matrix_2b' * q_sampled);
    
    % Extract parameters
    a1_minus_lambda1_2b = theta_lambda_2b(1);
    a2_minus_lambda2_2b = theta_lambda_2b(2);
    b0_2b = theta_lambda_2b(3);
    a1_2b = a1_minus_lambda1_2b + lambda1_2b;
    a2_2b = a2_minus_lambda2_2b + lambda2_2b;
    
    % Compute physical parameters
    c_est_2b = a1_2b / b0_2b;
    L_est_2b = g / a2_2b;
    m_est_2b = 1 / (b0_2b * L_est_2b^2);
    
    % Compute relative errors
    m_error_2b(i) = abs(m_est_2b - m_true) / m_true * 100;
    L_error_2b(i) = abs(L_est_2b - L_true) / L_true * 100;
    c_error_2b(i) = abs(c_est_2b - c_true) / c_true * 100;
    
    % Display progress
    fprintf('Processed input amplitude A0 = %.1f\n', A0);
end

%% Plot results
% Create separate plots for each parameter
figure('Position', [100, 100, 1200, 400]);

% Mass parameter error
subplot(1,3,1);
plot(input_amplitudes, m_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(input_amplitudes, m_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Input Amplitude A_0');
ylabel('Relative Error [%]');
title('Mass Parameter Error');
legend('Method 2a', 'Method 2b', 'Location', 'best');

% Length parameter error
subplot(1,3,2);
plot(input_amplitudes, L_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(input_amplitudes, L_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Input Amplitude A_0');
ylabel('Relative Error [%]');
title('Length Parameter Error');
legend('Method 2a', 'Method 2b', 'Location', 'best');

% Damping parameter error
subplot(1,3,3);
plot(input_amplitudes, c_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(input_amplitudes, c_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Input Amplitude A_0');
ylabel('Relative Error [%]');
title('Damping Parameter Error');
legend('Method 2a', 'Method 2b', 'Location', 'best');

sgtitle('Effect of Input Amplitude on Parameter Estimation Accuracy (T_s = 0.1 sec)');
