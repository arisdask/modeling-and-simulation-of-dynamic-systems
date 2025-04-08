% Modeling and Simulation of Dynamic Systems - Assignment 1
% Estimation of Unknown Parameters - Least Squares Method
% Author: Aristeidis Daskalopoulos - AEM:10640

%% Problem 3b: Effect of varying sampling period on parameter estimation accuracy
% Initialize Problem
clear; clc; close all;

addpath('src');
addpath('utils');

% True system parameters
m_true = 0.75; L_true = 1.25; c_true = 0.15;
g = 9.81; A0 = 4; omega = 2;

% Input torque function
inputTorque = @(t) A0 * sin(omega * t);

T_final = 20; % total simulation time [sec]
dt = 1e-4; % integration step

% Initial conditions
x0 = [0; 0]; % [q(0) q'(0)] = [0 0]

% Pendulum ODE function
pendulumODE = @(t, x, m, L, c, g, inputFunc) [
    x(2); % q' = dq/dt = x_1'
    (inputFunc(t) - c*x(2) - m*g*L*x(1))/(m*L^2) % q" = x_2'
];

% Generate "true" system data
[t_true, x_true] = solveODE(pendulumODE, 0:dt:T_final, x0, m_true, L_true, c_true, g, inputTorque);

% Define range of sampling periods to test
sampling_periods = 0.02:0.02:0.2;
num_periods = length(sampling_periods);

% Initialize arrays to store parameter estimation errors
m_error_2a = zeros(num_periods, 1);
L_error_2a = zeros(num_periods, 1);
c_error_2a = zeros(num_periods, 1);
m_error_2b = zeros(num_periods, 1);
L_error_2b = zeros(num_periods, 1);
c_error_2b = zeros(num_periods, 1);

% Define filter parameters
lambda1_2a = 0.2;     % For method 2a
lambda1_2b = 0.6;     % For method 2b
lambda2_2b = 0.08;    % For method 2b

% Loop through each sampling period
for i = 1:num_periods
    Ts = sampling_periods(i);
    
    % Sample the data at current Ts intervals
    sampled_indices = 1:round(Ts/dt):length(t_true);
    t_sampled = t_true(sampled_indices);
    x_sampled = x_true(sampled_indices, :);
    q_sampled = x_sampled(:, 1);
    q_dot_sampled = x_sampled(:, 2);
    u_sampled = inputTorque(t_sampled);
    
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
    fprintf('Processed sampling period Ts = %.3f s\n', Ts);
end

%% Plot results
% Create separate plots for each parameter
figure('Position', [100, 100, 1200, 400]);

% Mass parameter error
subplot(1,3,1);
plot(sampling_periods, m_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(sampling_periods, m_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Sampling Period T_s [sec]');
ylabel('Relative Error [%]');
title('Mass Parameter Error');
legend('Method 2a', 'Method 2b', 'Location', 'best');

% Length parameter error
subplot(1,3,2);
plot(sampling_periods, L_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(sampling_periods, L_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Sampling Period T_s [sec]');
ylabel('Relative Error [%]');
title('Length Parameter Error');
legend('Method 2a', 'Method 2b', 'Location', 'best');

% Damping parameter error
subplot(1,3,3);
plot(sampling_periods, c_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(sampling_periods, c_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Sampling Period T_s [sec]');
ylabel('Relative Error [%]');
title('Damping Parameter Error');
legend('Method 2a', 'Method 2b', 'Location', 'best');

sgtitle('Effect of Sampling Period on Parameter Estimation Accuracy');

% Create a combined plot
figure;
plot(sampling_periods, m_error_2a, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(sampling_periods, L_error_2a, 'b*-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(sampling_periods, c_error_2a, 'bd-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(sampling_periods, m_error_2b, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(sampling_periods, L_error_2b, 'r*-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(sampling_periods, c_error_2b, 'rd-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Sampling Period T_s [sec]');
ylabel('Relative Error [%]');
title('Combined Parameter Estimation Errors vs. Sampling Period');
legend('Mass (2a)', 'Length (2a)', 'Damping (2a)', 'Mass (2b)', 'Length (2b)', 'Damping (2b)', 'Location', 'best');
