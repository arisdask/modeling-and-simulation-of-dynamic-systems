%% Modeling and Simulation of Dynamic Systems - Assignment 2
% Author: Aristeidis Daskalopoulos - AEM:10640

% Mass-Spring-Damper System - Parameter Estimation using Lyapunov Method
% This script implements both Parallel and Mixed structure estimators for 
% the mass-spring-damper system.
%
% System:     m*d²x/dt² + b*dx/dt + k*x = u(t)
% True values: m = 1.315, b = 0.225, k = 0.725
% Input:      u(t) = 2.5*sin(t)

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

% Simulation parameters
T_final = 20;      % Simulation time (seconds)
tspan = [0 T_final];

%% Initial Conditions for Both Structures
% States: [x; dx/dt; x̂; dx̂/dt; θ1; θ2; θ3]
x0 = [0; 0];               % Initial true state [position; velocity]
x_hat0 = [0; 0];          % Initial estimated state
theta0 = [0.25; 0.5; 0.5];   % Initial parameter estimates

initial_conditions = [x0; x_hat0; theta0];

%% Parallel Structure Simulation
[t_parallel, y_parallel] = ode45( ...
    @(t,y) parallelStructureDynamics(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3), ...
    tspan, ...
    initial_conditions);

% Extract results
x_parallel = y_parallel(:, 1);           % True position
xdot_parallel = y_parallel(:, 2);        % True velocity
x_hat_parallel = y_parallel(:, 3);       % Estimated position
xdot_hat_parallel = y_parallel(:, 4);    % Estimated velocity
theta_parallel = y_parallel(:, 5:7);     % Parameter estimates

% Compute physical parameters
m_est_parallel = 1 ./ theta_parallel(:, 3);
b_est_parallel = theta_parallel(:, 1) .* m_est_parallel;
k_est_parallel = theta_parallel(:, 2) .* m_est_parallel;

% Compute estimation errors
e_x_parallel = x_parallel - x_hat_parallel;
e_xdot_parallel = xdot_parallel - xdot_hat_parallel;

%% Mixed Structure Simulation
[t_mixed, y_mixed] = ode45( ...
    @(t,y) mixedStructureDynamics(t, y, m_true, b_true, k_true, gamma1, gamma2, gamma3, theta_m), ...
    tspan, ...
    initial_conditions);

% Extract results
x_mixed = y_mixed(:, 1);           % True position
xdot_mixed = y_mixed(:, 2);        % True velocity
x_hat_mixed = y_mixed(:, 3);       % Estimated position
xdot_hat_mixed = y_mixed(:, 4);    % Estimated velocity
theta_mixed = y_mixed(:, 5:7);     % Parameter estimates

% Compute physical parameters
m_est_mixed = 1 ./ theta_mixed(:, 3);
b_est_mixed = theta_mixed(:, 1) .* m_est_mixed;
k_est_mixed = theta_mixed(:, 2) .* m_est_mixed;

% Compute estimation errors
e_x_mixed = x_mixed - x_hat_mixed;
e_xdot_mixed = xdot_mixed - xdot_hat_mixed;

%% Plotting Results

% Figure 1: Parallel Structure Results
figure('Name', 'Parallel Structure Results');

% States and Error
subplot(2,1,1);
plot(t_parallel, x_parallel, 'b-', t_parallel, x_hat_parallel, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Position x(t) [m]');
title('Position: x(t) and x̂(t)');
legend('True x(t)', 'Estimated x̂(t)');

subplot(2,1,2);
plot(t_parallel, e_x_parallel, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error e_x(t) [m]');
title('Position Error: e_x(t) = x(t) - x̂(t)');

% Parameter Estimates
figure('Name', 'Parallel Structure Results');
plot(t_parallel, m_est_parallel, 'b-', ...
     t_parallel, b_est_parallel, 'r--', ...
     t_parallel, k_est_parallel, 'g:', ...
     [0 T_final], [m_true m_true], 'b--', ...
     [0 T_final], [b_true b_true], 'r:', ...
     [0 T_final], [k_true k_true], 'g-.', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Parameter Values');
title('Parameter Estimates');
legend('m̂(t)', 'b̂(t)', 'k̂(t)', 'm_{true}', 'b_{true}', 'k_{true}');


sgtitle('Parallel Structure - Parameter Estimation Results', 'FontSize', 14);

% Figure 2: Mixed Structure Results
figure('Name', 'Mixed Structure Results');

% States and Error
subplot(2,1,1);
plot(t_mixed, x_mixed, 'b-', t_mixed, x_hat_mixed, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Position x(t) [m]');
title('Position: x(t) and x̂(t)');
legend('True x(t)', 'Estimated x̂(t)');

subplot(2,1,2);
plot(t_mixed, e_x_mixed, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error e_x(t) [m]');
title('Position Error: e_x(t) = x(t) - x̂(t)');

% Parameter Estimates
figure('Name', 'Mixed Structure Results');
plot(t_mixed, m_est_mixed, 'b-', ...
     t_mixed, b_est_mixed, 'r--', ...
     t_mixed, k_est_mixed, 'g:', ...
     [0 T_final], [m_true m_true], 'b--', ...
     [0 T_final], [b_true b_true], 'r:', ...
     [0 T_final], [k_true k_true], 'g-.', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Parameter Values');
title('Parameter Estimates');
legend('m̂(t)', 'b̂(t)', 'k̂(t)', 'm_{true}', 'b_{true}', 'k_{true}');

%% Display Final Results
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
