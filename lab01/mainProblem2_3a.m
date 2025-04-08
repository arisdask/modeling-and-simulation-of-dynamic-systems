% Modeling and Simulation of Dynamic Systems - Assignment 1
% Estimation of Unknown Parameters - Least Squares Method
% Author: Aristeidis Daskalopoulos - AEM:10640

%% Initialize Problem
clear; clc; close all;

addpath('src');
addpath('utils');

% True system parameters
m_true = 0.75; L_true = 1.25; c_true = 0.15;
g = 9.81; A0 = 4; omega = 2;

% Input Torque function
inputTorque = @(t) A0 * sin(omega * t);

T_final = 20;   % total simulation time [sec]
dt      = 1e-4; % integration step
Ts      = 0.1;  % sampling period for parameter estimation [sec]

time_continuous = 0:dt:T_final;
time_sampled = 0:Ts:T_final;

% Initial conditions
x0 = [0; 0];  % [q(0) q'(0)] = [0 0]

% Pendulum ODE function (state x = [x_1 x_2]^T)
pendulumODE = @(t, x, m, L, c, g, inputFunc) [
    x(2);                                           % q' = dq/dt = x_1'
    (inputFunc(t) - c*x(2) - m*g*L*x(1))/(m*L^2)    % q" = x_2'
];

%% Generate true/"real" system data (Problem 1)
[t_true, x_true] = solveODE(pendulumODE, time_continuous, x0, m_true, L_true, c_true, g, inputTorque);

% Sample the data at Ts intervals
sampled_indices = 1:Ts/dt:length(t_true);

t_sampled       = t_true(sampled_indices);
x_sampled       = x_true(sampled_indices, :);
q_sampled       = x_sampled(:, 1);  % x_1 = q
q_dot_sampled   = x_sampled(:, 2);  % x_2 = q'
u_sampled       = inputTorque(t_sampled);

%% Problem 2a: Estimate parameters when q(t), q'(t), and u(t) are measurable
% Define filter parameter lambda1 (pole of the first-order filter)
lambda1 = 0.2;  % The filter is defined as 1/Lambda(s) where Lambda(s) = s + lambda1

% Create the filtered signals according to the mathematical analysis
% The filter transfer functions are:
% ζ1 = -1/Lambda(s) * q_dot
% ζ2 = -1/Lambda(s) * q
% ζ3 = +1/Lambda(s) * u

% Create the filter transfer functions
s = tf('s');
Lambda_s = s + lambda1;

% Create the filtered signals
zeta1 = lsim(-1/Lambda_s, q_dot_sampled, t_sampled);   % -1/Lambda(s) * q_dot
zeta2 = lsim(-1/Lambda_s, q_sampled, t_sampled);       % -1/Lambda(s) * q
zeta3 = lsim(1/Lambda_s, u_sampled, t_sampled);        % +1/Lambda(s) * u

% From the mathematical analysis, we have q_dot = θ_λ^T * ζ where θ_λ = [a1-λ1, a2, b0]^T
% Create the ζ matrix for least squares regression
zeta_matrix = [zeta1, zeta2, zeta3];

% Apply least squares to solve for θ_λ
theta_lambda = (zeta_matrix' * zeta_matrix) \ (zeta_matrix' * q_dot_sampled);

% Extract the parameters
a1_minus_lambda1 = theta_lambda(1);
a2 = theta_lambda(2);
b0 = theta_lambda(3);

% Calculate a1
a1 = a1_minus_lambda1 + lambda1;

% Extract the physical parameters m, L, c
c_est_2a          = a1 / b0;                % c = a1 * mL^2 and b0 = 1 / mL^2
L_est_2a          = g / a2;                 % a2 = g / L
m_est_2a          = 1 / (b0 * L_est_2a^2);  % b0 = 1/(m*L^2)

fprintf('\nProblem 2a - Estimated Parameters:\n');
printEstimationResults('m', m_est_2a, m_true);
printEstimationResults('L', L_est_2a, L_true);
printEstimationResults('c', c_est_2a, c_true);

% Simulate and compute errors
[t_est_2a, x_est_2a, e_q_2a, e_q_dot_2a] = ...
    simulateAndComputeErrors( ...
    pendulumODE, time_continuous, x0, m_est_2a, L_est_2a, c_est_2a, g, inputTorque, t_true, x_true ...
    );

% Plot results
plotResults(t_true, x_true, t_est_2a, x_est_2a, e_q_2a, e_q_dot_2a, ...
    'Problem 2a: Parameter Estimation');

% Print RMS errors
printRMSErrors(e_q_2a, e_q_dot_2a);

%% Problem 2b: Estimate parameters when only q(t) and u(t) are measurable
% Define filter parameters lambda1 and lambda2
% These values should be chosen to place the poles in the left half-plane
% The filter is defined as 1/Lambda(s) where Lambda(s) = s^2 + lambda1*s + lambda2
lambda1 = 0.6;
lambda2 = 0.08;

% Create the filtered signals according to the mathematical analysis
% The filter transfer functions are:
% ζ1 = -s/Lambda(s) * q
% ζ2 = -1/Lambda(s) * q
% ζ3 = +1/Lambda(s) * u

% Create the filter transfer functions
s = tf('s');
Lambda_s = s^2 + lambda1*s + lambda2;

% Create the filtered signals
zeta1 = lsim(-s/Lambda_s, q_sampled, t_sampled);   % -s/Lambda(s) * q
zeta2 = lsim(-1/Lambda_s, q_sampled, t_sampled);   % -1/Lambda(s) * q
zeta3 = lsim(1/Lambda_s, u_sampled, t_sampled);    % +1/Lambda(s) * u

% From the mathematical analysis, we have q = θ_λ^T * ζ where θ_λ = [a1-λ1, a2-λ2, b0]^T
% Create the ζ matrix for least squares regression
zeta_matrix = [zeta1, zeta2, zeta3];

% Apply least squares to solve for θ_λ
theta_lambda = (zeta_matrix' * zeta_matrix) \ (zeta_matrix' * q_sampled);

% Extract the original parameters
a1_minus_lambda1 = theta_lambda(1);
a2_minus_lambda2 = theta_lambda(2);
b0 = theta_lambda(3);

% Calculate the actual parameters a1 and a2
a1 = a1_minus_lambda1 + lambda1;
a2 = a2_minus_lambda2 + lambda2;

% Extract the physical parameters m, L, c (once again :/)
c_est_2b          = a1 / b0;                % c = a1 * mL^2 and b0 = 1 / mL^2
L_est_2b          = g / a2;                 % a2 = g / L
m_est_2b          = 1 / (b0 * L_est_2b^2);  % b0 = 1/(m*L^2)

fprintf('\nProblem 2b - Estimated Parameters:\n');
printEstimationResults('m', m_est_2b, m_true);
printEstimationResults('L', L_est_2b, L_true);
printEstimationResults('c', c_est_2b, c_true);

% Simulate and compute errors
[t_est_2b, x_est_2b, e_q_2b, e_q_dot_2b] = ...
    simulateAndComputeErrors( ...
    pendulumODE, time_continuous, x0, m_est_2b, L_est_2b, c_est_2b, g, inputTorque, t_true, x_true ...
    );

% Plot results
plotResults(t_true, x_true, t_est_2b, x_est_2b, e_q_2b, e_q_dot_2b, ...
    'Problem 2b: Parameter Estimation');

% Print RMS errors
printRMSErrors(e_q_2b, e_q_dot_2b);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Problem 3a: Parameter estimation with noisy measurements
% Noise Levels
noise_level_q = 0.05;      % Standard deviation for angle noise
noise_level_q_dot = 0.05;  % Standard deviation for angular velocity noise

% Add noise to sampled data
[q_noisy, q_dot_noisy, u_sampled] = addNoise(q_sampled, q_dot_sampled, u_sampled, noise_level_q, noise_level_q_dot);

% Plot original vs noisy data for comparison
figure;
subplot(2,1,1);
plot(t_sampled, q_sampled, 'b-', t_sampled, q_noisy, 'r--', 'LineWidth', 1.5);
title('Angle q(t) with and without noise');
xlabel('Time [sec]');
ylabel('Angle [rad]');
legend('Clean', 'Noisy');
grid on;

subplot(2,1,2);
plot(t_sampled, q_dot_sampled, 'b-', t_sampled, q_dot_noisy, 'r--', 'LineWidth', 1.5);
title('Angular Velocity q̇(t) with and without noise');
xlabel('Time [sec]');
ylabel('Angular Velocity [rad/s]');
legend('Clean', 'Noisy');
grid on;

%% Problem 3a.1: Estimate parameters with noisy data using approach from Problem 2a
% Define filter parameter lambda1 (pole of the first-order filter)
lambda1 = 0.2;

% Create the filter transfer functions
s = tf('s');
Lambda_s = s + lambda1;

% Create the filtered signals with noisy data
zeta1_noisy = lsim(-1/Lambda_s, q_dot_noisy, t_sampled);   % -1/Lambda(s) * q_dot
zeta2_noisy = lsim(-1/Lambda_s, q_noisy, t_sampled);       % -1/Lambda(s) * q
zeta3_noisy = lsim(1/Lambda_s, u_sampled, t_sampled);      % 1/Lambda(s) * u

% Create the ζ matrix for least squares regression
zeta_matrix_noisy = [zeta1_noisy, zeta2_noisy, zeta3_noisy];

% Apply least squares to solve for θ_λ
theta_lambda_noisy = (zeta_matrix_noisy' * zeta_matrix_noisy) \ (zeta_matrix_noisy' * q_dot_noisy);

% Extract the parameters
a1_minus_lambda1_noisy = theta_lambda_noisy(1);
a2_noisy = theta_lambda_noisy(2);
b0_noisy = theta_lambda_noisy(3);

% Calculate a1
a1_noisy = a1_minus_lambda1_noisy + lambda1;

% From these, extract the physical parameters m, L, c
c_est_noisy_3a1 = a1_noisy / b0_noisy;
L_est_noisy_3a1 = g / a2_noisy;
m_est_noisy_3a1 = 1 / (b0_noisy * L_est_noisy_3a1^2);

fprintf('\nProblem 3a.1 - Estimated Parameters with Noisy Data (Method 2a):\n');
compareNoiseImpact('m', m_est_2a, m_est_noisy_3a1, m_true);
compareNoiseImpact('L', L_est_2a, L_est_noisy_3a1, L_true);
compareNoiseImpact('c', c_est_2a, c_est_noisy_3a1, c_true);

% Simulate the system with estimated parameters from noisy data
[t_est_noisy_3a1, x_est_noisy_3a1, e_q_noisy_3a1, ~] = simulateAndComputeErrors(pendulumODE, time_continuous, x0, m_est_noisy_3a1, L_est_noisy_3a1, c_est_noisy_3a1, g, inputTorque, t_true, x_true);

% Plot angle results (q(t)) with noisy estimates vs true data
figure;
subplot(2,1,1);
plot(t_true, x_true(:,1), 'b-', t_est_noisy_3a1, x_est_noisy_3a1(:,1), 'r--', 'LineWidth', 1.5);
title('Angle q(t): True vs. Estimated from Noisy Data');
xlabel('Time [sec]');
ylabel('Angle [rad]');
legend('True', 'Estimated from Noisy Data');
grid on;

subplot(2,1,2);
plot(t_est_noisy_3a1, e_q_noisy_3a1, 'k-', 'LineWidth', 1.5);
title('Error in Angle e_q(t) with Noisy Estimation');
xlabel('Time [sec]');
ylabel('Error [rad]');
grid on;

sgtitle('Problem 3a.1: Parameter Estimation from Noisy Data (Method 2a)');

% Calculate RMS errors for angle only
rms_e_q_noisy_3a1 = sqrt(mean(e_q_noisy_3a1.^2));
fprintf('\nRMS Error with Noisy Data (Method 2a):\n');
fprintf('RMS e_q = %.6f\n', rms_e_q_noisy_3a1);

%% Problem 3a.2: Estimate parameters with noisy data using approach from Problem 2b
% Define filter parameters lambda1 and lambda2
lambda1 = 0.6;
lambda2 = 0.08;

% Create the filter transfer functions
s = tf('s');
Lambda_s = s^2 + lambda1*s + lambda2;

% Create the filtered signals with noisy data
zeta1_noisy = lsim(-s/Lambda_s, q_noisy, t_sampled);   % -s/Lambda(s) * q
zeta2_noisy = lsim(-1/Lambda_s, q_noisy, t_sampled);   % -1/Lambda(s) * q
zeta3_noisy = lsim(1/Lambda_s, u_sampled, t_sampled);  % 1/Lambda(s) * u

% Create the ζ matrix for least squares regression
zeta_matrix_noisy = [zeta1_noisy, zeta2_noisy, zeta3_noisy];

% Apply least squares to solve for θ_λ
theta_lambda_noisy = (zeta_matrix_noisy' * zeta_matrix_noisy) \ (zeta_matrix_noisy' * q_noisy);

% Extract the original parameters
a1_minus_lambda1_noisy = theta_lambda_noisy(1);
a2_minus_lambda2_noisy = theta_lambda_noisy(2);
b0_noisy = theta_lambda_noisy(3);

% Calculate the actual parameters a1 and a2
a1_noisy = a1_minus_lambda1_noisy + lambda1;
a2_noisy = a2_minus_lambda2_noisy + lambda2;

% From these, extract the physical parameters m, L, c
c_est_noisy_3a2 = a1_noisy / b0_noisy;
L_est_noisy_3a2 = g / a2_noisy;
m_est_noisy_3a2 = 1 / (b0_noisy * L_est_noisy_3a2^2);

fprintf('\nProblem 3a.2 - Estimated Parameters with Noisy Data (Method 2b):\n');
compareNoiseImpact('m', m_est_2b, m_est_noisy_3a2, m_true);
compareNoiseImpact('L', L_est_2b, L_est_noisy_3a2, L_true);
compareNoiseImpact('c', c_est_2b, c_est_noisy_3a2, c_true);

% Simulate the system with estimated parameters from noisy data
[t_est_noisy_3a2, x_est_noisy_3a2, e_q_noisy_3a2, ~] = simulateAndComputeErrors(pendulumODE, time_continuous, x0, m_est_noisy_3a2, L_est_noisy_3a2, c_est_noisy_3a2, g, inputTorque, t_true, x_true);

% Plot angle results (q(t)) with noisy estimates vs true data
figure;
subplot(2,1,1);
plot(t_true, x_true(:,1), 'b-', t_est_noisy_3a2, x_est_noisy_3a2(:,1), 'r--', 'LineWidth', 1.5);
title('Angle q(t): True vs. Estimated from Noisy Data');
xlabel('Time [sec]');
ylabel('Angle [rad]');
legend('True', 'Estimated from Noisy Data');
grid on;

subplot(2,1,2);
plot(t_est_noisy_3a2, e_q_noisy_3a2, 'k-', 'LineWidth', 1.5);
title('Error in Angle e_q(t) with Noisy Estimation');
xlabel('Time [sec]');
ylabel('Error [rad]');
grid on;

sgtitle('Problem 3a.2: Parameter Estimation from Noisy Data (Method 2b)');

% Calculate RMS errors for angle only
rms_e_q_noisy_3a2 = sqrt(mean(e_q_noisy_3a2.^2));
fprintf('\nRMS Error with Noisy Data (Method 2b):\n');
fprintf('RMS e_q = %.6f\n', rms_e_q_noisy_3a2);

%% Compare all methods side by side
figure;
plot(t_true, x_true(:,1), 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_est_2a, x_est_2a(:,1), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Clean 2a');
plot(t_est_2b, x_est_2b(:,1), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Clean 2b');
plot(t_est_noisy_3a1, x_est_noisy_3a1(:,1), 'r:', 'LineWidth', 1.5, 'DisplayName', 'Noisy 2a');
plot(t_est_noisy_3a2, x_est_noisy_3a2(:,1), 'm:', 'LineWidth', 1.5, 'DisplayName', 'Noisy 2b');
title('Comparison of All Estimation Methods - Angle q(t)');
xlabel('Time [sec]');
ylabel('Angle [rad]');
legend('Location', 'best');
grid on;

% Print comparison table of parameter estimates
fprintf('\nComparison of Parameter Estimates:\n');
fprintf('Parameter | True     | Clean 2a  | Clean 2b  | Noisy 2a  | Noisy 2b\n');
fprintf('------------------------------------------------------------------\n');
fprintf('m         | %.4f    | %.4f    | %.4f    | %.4f    | %.4f\n', ...
    m_true, m_est_2a, m_est_2b, m_est_noisy_3a1, m_est_noisy_3a2);
fprintf('L         | %.4f    | %.4f    | %.4f    | %.4f    | %.4f\n', ...
    L_true, L_est_2a, L_est_2b, L_est_noisy_3a1, L_est_noisy_3a2);
fprintf('c         | %.4f    | %.4f    | %.4f    | %.4f    | %.4f\n', ...
    c_true, c_est_2a, c_est_2b, c_est_noisy_3a1, c_est_noisy_3a2);
