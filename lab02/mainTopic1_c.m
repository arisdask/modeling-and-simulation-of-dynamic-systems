%% Modeling and Simulation of Dynamic Systems - Assignment 2
% Author: Aristeidis Daskalopoulos - AEM:10640

% Mass-Spring-Damper System - Parameter Estimation using Lyapunov Method
% with Noise
% This script implements both Parallel and Mixed structure estimators for 
% the mass-spring-damper system.
%
% System:     m*d²x/dt² + b*dx/dt + k*x = u(t)
% True values: m = 1.315, b = 0.225, k = 0.725
% Input:      u(t) = 2.5*sin(t)
% Noise:      n(t) = n0*sin(2πf0*t), f0 = 20Hz

clear; clc; close all;

addpath('src');
addpath('utils');

%% System Parameters
m_true = 1.315;    % True mass value
b_true = 0.225;    % True damping coefficient
k_true = 0.725;    % True spring stiffness

% Learning rates for both structures
gamma1 = 2.8;      % For θ1 (related to b/m)
gamma2 = 0.06;     % For θ2 (related to k/m)
gamma3 = 2.0;      % For θ3 (related to 1/m)

% Mixed structure specific parameter
theta_m = 60;     % Matrix C coefficient

% Noise parameters
f0 = 20;          % Noise frequency [Hz]
n0_base = 0.25;   % Base noise amplitude

% Simulation parameters
T_final = 20;      % Simulation time (seconds)
tspan = [0 T_final];

%% Initial Conditions for Both Structures
% States: [x; dx/dt; x̂; dx̂/dt; θ1; θ2; θ3]
x0 = [0; 0];               % Initial true state [position; velocity]
x_hat0 = [0; 0];          % Initial estimated state
theta0 = [0.25; 0.5; 0.5];   % Initial parameter estimates

initial_conditions = [x0; x_hat0; theta0];

%% Part 1: Simulation with base noise amplitude (n0 = 0.25)
% Parallel Structure Simulation with noise
[t_parallel_noise, y_parallel_noise] = ode45( ...
    @(t,y) parallelStructureDynamicsNoise(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3, n0_base, f0), ...
    tspan, ...
    initial_conditions);

% Mixed Structure Simulation with noise
[t_mixed_noise, y_mixed_noise] = ode45( ...
    @(t,y) mixedStructureDynamicsNoise(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3, theta_m, n0_base, f0), ...
    tspan, ...
    initial_conditions);

% Extract and process results for both structures
[x_parallel, xdot_parallel, x_hat_parallel, xdot_hat_parallel, ...
    m_est_parallel, b_est_parallel, k_est_parallel, e_x_parallel, x_measured_parallel] = ...
    extractResultsWithNoise(y_parallel_noise, t_parallel_noise, n0_base, f0);

[x_mixed, xdot_mixed, x_hat_mixed, xdot_hat_mixed, ...
    m_est_mixed, b_est_mixed, k_est_mixed, e_x_mixed, x_measured_mixed] = ...
    extractResultsWithNoise(y_mixed_noise, t_mixed_noise, n0_base, f0);

% Plot results for base noise case
plotResultsWithNoise(t_parallel_noise, x_parallel, x_hat_parallel, x_measured_parallel, ...
    e_x_parallel, m_est_parallel, b_est_parallel, k_est_parallel, ...
    m_true, b_true, k_true, T_final, ...
    'Parallel Structure');

plotResultsWithNoise(t_mixed_noise, x_mixed, x_hat_mixed, x_measured_mixed, ...
    e_x_mixed, m_est_mixed, b_est_mixed, k_est_mixed, ...
    m_true, b_true, k_true, T_final, ...
    'Mixed Structure');

%% Part 2: Study effect of varying noise amplitude
% Define range of noise amplitudes to test
n0_values = 0.1:0.1:1;
num_n0 = length(n0_values);

% Initialize arrays to store parameter errors
m_error_parallel = zeros(num_n0, 1);
b_error_parallel = zeros(num_n0, 1);
k_error_parallel = zeros(num_n0, 1);
m_error_mixed = zeros(num_n0, 1);
b_error_mixed = zeros(num_n0, 1);
k_error_mixed = zeros(num_n0, 1);

% Run simulations for each noise amplitude
for i = 1:num_n0
    n0 = n0_values(i);
    
    % Parallel Structure
    [t_p, y_p] = ode45(@(t,y) parallelStructureDynamicsNoise(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3, n0, f0), ...
        tspan, initial_conditions);
    [~, ~, ~, ~, m_est_p, b_est_p, k_est_p, ~, ~] = extractResultsWithNoise(y_p, t_p, n0, f0);
    
    % Mixed Structure
    [t_m, y_m] = ode45(@(t,y) mixedStructureDynamicsNoise(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3, theta_m, n0, f0), ...
        tspan, initial_conditions);
    [~, ~, ~, ~, m_est_m, b_est_m, k_est_m, ~, ~] = extractResultsWithNoise(y_m, t_m, n0, f0);
    
    % Compute relative errors (using last 5 points average)
    m_error_parallel(i) = 100*abs(mean(m_est_p(end-5:end))-m_true)/m_true;
    b_error_parallel(i) = 100*abs(mean(b_est_p(end-5:end))-b_true)/b_true;
    k_error_parallel(i) = 100*abs(mean(k_est_p(end-5:end))-k_true)/k_true;
    
    m_error_mixed(i) = 100*abs(mean(m_est_m(end-5:end))-m_true)/m_true;
    b_error_mixed(i) = 100*abs(mean(b_est_m(end-5:end))-b_true)/b_true;
    k_error_mixed(i) = 100*abs(mean(k_est_m(end-5:end))-k_true)/k_true;
end

% Plot parameter errors vs noise amplitude in separate figures
% Mass parameter error
figure('Name', 'Mass Parameter Error vs Noise Amplitude');
subplot(2,1,1);
plot(n0_values, m_error_parallel, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Noise Amplitude n_0');
ylabel('Mass Error [%]');
title('Parallel Structure - Mass Parameter Error');

subplot(2,1,2);
plot(n0_values, m_error_mixed, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Noise Amplitude n_0');
ylabel('Mass Error [%]');
title('Mixed Structure - Mass Parameter Error');
sgtitle('Mass Parameter Estimation Error vs Noise Amplitude');

% Damping parameter error
figure('Name', 'Damping Parameter Error vs Noise Amplitude');
subplot(2,1,1);
plot(n0_values, b_error_parallel, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Noise Amplitude n_0');
ylabel('Damping Error [%]');
title('Parallel Structure - Damping Parameter Error');

subplot(2,1,2);
plot(n0_values, b_error_mixed, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Noise Amplitude n_0');
ylabel('Damping Error [%]');
title('Mixed Structure - Damping Parameter Error');
sgtitle('Damping Parameter Estimation Error vs Noise Amplitude');

% Spring constant parameter error
figure('Name', 'Spring Constant Parameter Error vs Noise Amplitude');
subplot(2,1,1);
plot(n0_values, k_error_parallel, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Noise Amplitude n_0');
ylabel('Spring Constant Error [%]');
title('Parallel Structure - Spring Constant Parameter Error');

subplot(2,1,2);
plot(n0_values, k_error_mixed, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Noise Amplitude n_0');
ylabel('Spring Constant Error [%]');
title('Mixed Structure - Spring Constant Parameter Error');
sgtitle('Spring Constant Parameter Estimation Error vs Noise Amplitude');

%% Display Final Results for base noise case
fprintf('\nResults with noise amplitude n0 = %.3f:\n', n0_base);
fprintf('\nParallel Structure Results:\n');
fprintf('True parameters:  m = %.3f, b = %.3f, k = %.3f\n', m_true, b_true, k_true);
fprintf('Final estimates: m = %.3f, b = %.3f, k = %.3f\n', ...
    mean(m_est_parallel(end-5:end)), ...
    mean(b_est_parallel(end-5:end)), ...
    mean(k_est_parallel(end-5:end)));
fprintf('Relative errors: m = %.2f%%, b = %.2f%%, k = %.2f%%\n', ...
    100*abs(mean(m_est_parallel(end-5:end))-m_true)/m_true, ...
    100*abs(mean(b_est_parallel(end-5:end))-b_true)/b_true, ...
    100*abs(mean(k_est_parallel(end-5:end))-k_true)/k_true);

fprintf('\nMixed Structure Results:\n');
fprintf('True parameters:  m = %.3f, b = %.3f, k = %.3f\n', m_true, b_true, k_true);
fprintf('Final estimates: m = %.3f, b = %.3f, k = %.3f\n', ...
    mean(m_est_mixed(end-5:end)), ...
    mean(b_est_mixed(end-5:end)), ...
    mean(k_est_mixed(end-5:end)));
fprintf('Relative errors: m = %.2f%%, b = %.2f%%, k = %.2f%%\n', ...
    100*abs(mean(m_est_mixed(end-5:end))-m_true)/m_true, ...
    100*abs(mean(b_est_mixed(end-5:end))-b_true)/b_true, ...
    100*abs(mean(k_est_mixed(end-5:end))-k_true)/k_true);
