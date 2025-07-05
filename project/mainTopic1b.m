%% Modeling and Simulation of Dynamic Systems - Project - Topic 1b
% Author: Aristeidis Daskalopoulos - AEM:10640
% Implementation of gradient method with σ-modification for biased systems

clear; clc; close all;
addpath('src\');
addpath('utils\');

%% True System Parameters
A_true = [-2.15,  0.25;
          -0.75, -2.0];

B_true = [0; 1.5];

%% Bias Error Parameters
% Modeling the bias as bounded sinusoidal disturbance
omega_bar = 0.3;  % Maximum bound for ||ω(t)||
omega_freq = [0.7, 1.2];  % Frequencies for bias components

%% Gradient Method Parameters

% Filter Parameters
lambda = 2.0;  % Filter pole (must be > 0 for stability)

% Learning rate matrix (diagonal)
gamma_vals = [1; 1; 1; 1; 1; 1];  % Individual learning rates
Gamma = diag(gamma_vals);

% σ-modification parameters
M = 3.4;        % Bound for parameter estimates
sigma = 1;   % σ-modification coefficient

%% Initial Conditions
% True system initial state
x0 = [0; 0];

% Initial parameter estimates (feasible and within bounds)
% [a11+λ, a12, a21, a22+λ, b1, b2]
theta_lambda0 = [-0.5; 3; 0; 0; 5; 5];

% Initialize filter states (all start at zero)
filter_states0 = zeros(3, 1);  % 3 filters (for x1, x2 and u)

% Combined initial conditions: [x; theta_lambda; filter_states]
initial_conditions = [x0; theta_lambda0; filter_states0];

%% Run Simulation

% Simulation Parameters
T_final = 30;  % Simulation time (seconds)
tspan = [0 T_final];

[t, x_sim] = ode45(@(t, x) ...
    biasedSystemDynamics(t, x, A_true, B_true, lambda, Gamma, M, sigma, omega_bar, omega_freq), ...
    tspan, initial_conditions);

%% Extract Results
n_states = 2;
n_params = 6;
n_filters = 12;

% Extract trajectories
x_true            =  x_sim(:, 1:n_states);
theta_lambda_est  =  x_sim(:, (n_states+1):(n_states+n_params));
filter_states     =  x_sim(:, (n_states+n_params+1):end);

% Convert θ_λ back to original parameters
theta_est = theta_lambda_est;
theta_est(:, 1) = theta_lambda_est(:, 1) - lambda;  % a11 = (a11 + λ) - λ
theta_est(:, 4) = theta_lambda_est(:, 4) - lambda;  % a22 = (a22 + λ) - λ

% Reconstruct estimated matrices
A_est = zeros(2, 2, length(t));
B_est = zeros(2, 1, length(t));

for i = 1:length(t)
    A_est(:, :, i) = [theta_est(i, 1), theta_est(i, 2);
                      theta_est(i, 3), theta_est(i, 4)];
    B_est(:, :, i) = [theta_est(i, 5); theta_est(i, 6)];
end

% Compute estimated states using estimated parameters
x_hat = zeros(size(x_true));
for i = 1:length(t)
    % Input signal
    u = 2.5 * sin(t(i)) + 1.5 * cos(1.5*t(i)) + 0.8 * sin(3*t(i));
    
    % Compute regressor matrix Φ(t)
    % Extract filtered signals from filter states
    x1_filt = filter_states(i, 1);
    x2_filt = filter_states(i, 2);
    u_filt  = filter_states(i, 3);
    
    Phi = [x1_filt, x2_filt, 0,       0,       u_filt, 0;
           0,       0,       x1_filt, x2_filt, 0,      u_filt];
    
    % Estimated output: ŷ = Φ * θ_λ
    x_hat(i, :) = (Phi * theta_lambda_est(i, :)')';
end

% Compute estimation error
e_x = x_true - x_hat;

%% Plot Results
% State trajectories and estimation
figure('Name', 'State Trajectories', 'NumberTitle', 'off');
subplot(2,2,1);
plot(t, x_true(:,1), 'b-', t, x_hat(:,1), 'r--', 'LineWidth', 1.5);
title('State x_1(t)');
xlabel('Time [s]');
ylabel('x_1');
legend('True', 'Estimated', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t, x_true(:,2), 'b-', t, x_hat(:,2), 'r--', 'LineWidth', 1.5);
title('State x_2(t)');
xlabel('Time [s]');
ylabel('x_2');
legend('True', 'Estimated', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(t, e_x(:,1), 'g-', 'LineWidth', 1.5);
title('Estimation Error e_1(t)');
xlabel('Time [s]');
ylabel('$e_1 = x_1 - \hat{x}_1$', 'Interpreter', 'latex');
grid on;

subplot(2,2,4);
plot(t, e_x(:,2), 'g-', 'LineWidth', 1.5);
title('Estimation Error e_2(t)');
xlabel('Time [s]');
ylabel('$e_2 = x_2 - \hat{x}_2$', 'Interpreter', 'latex');
grid on;

% Estimation error norm
figure('Name', 'Estimation Error Norm', 'NumberTitle', 'off');
plot(t, sqrt(e_x(:,1).^2 + e_x(:,2).^2), 'm-', 'LineWidth', 1.5);
title('Estimation Error Norm ||e_x||');
xlabel('Time [s]');
ylabel('||e_x||');
grid on;

% Parameter estimates for matrix A
figure('Name', 'Matrix A Parameter Estimates', 'NumberTitle', 'off');
subplot(2,2,1);
plot(t, theta_est(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(1,1)*ones(size(t)), 'b--', 'LineWidth', 1.5);
plot(t, -3*ones(size(t)), 'k:', 'LineWidth', 1);
plot(t, -1*ones(size(t)), 'k:', 'LineWidth', 1);
title('Parameter a_{11}');
xlabel('Time [s]');
ylabel('a_{11}');
legend('Estimated', 'True', 'Constraints', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t, theta_est(:,2), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(1,2)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter a_{12}');
xlabel('Time [s]');
ylabel('a_{12}');
legend('Estimated', 'True', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(t, theta_est(:,3), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(2,1)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter a_{21}');
xlabel('Time [s]');
ylabel('a_{21}');
legend('Estimated', 'True', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(t, theta_est(:,4), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(2,2)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter a_{22}');
xlabel('Time [s]');
ylabel('a_{22}');
legend('Estimated', 'True', 'Location', 'best');
grid on;

% Parameter estimates for matrix B
figure('Name', 'Matrix B Parameter Estimates', 'NumberTitle', 'off');
subplot(1,2,1);
plot(t, theta_est(:,5), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, B_true(1)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter b_1');
xlabel('Time [s]');
ylabel('b_1');
legend('Estimated', 'True', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(t, theta_est(:,6), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, B_true(2)*ones(size(t)), 'b--', 'LineWidth', 1.5);
plot(t, ones(size(t)), 'k:', 'LineWidth', 1);
title('Parameter b_2');
xlabel('Time [s]');
ylabel('b_2');
legend('Estimated', 'True', 'Constraint b_2 \geq 1', 'Location', 'best');
grid on;

% Parameter estimate norms for σ-modification visualization
figure('Name', 'Parameter Norms and σ-Modification', 'NumberTitle', 'off');
theta_norms = sqrt(sum(theta_lambda_est.^2, 2));
plot(t, theta_norms, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, M*ones(size(t)), 'r--', 'LineWidth', 1.5);
plot(t, 2*M*ones(size(t)), 'r:', 'LineWidth', 1.5);
title('Parameter Estimate Norm ||θ_λ||');
xlabel('Time [s]');
ylabel('||θ_λ||');
legend('||θ_λ||', 'M', '2M', 'Location', 'best');
grid on;

%% Display final estimates
fprintf('=== PART B RESULTS ===\n');
fprintf('Final Parameter Estimates:\n');
fprintf('a11: %.4f (True: %.4f - Error: %.2f%%)\n', ...
    theta_est(end, 1), A_true(1,1), 100 * abs(theta_est(end, 1) - A_true(1,1))/abs(A_true(1,1)));
fprintf('a12: %.4f (True: %.4f - Error: %.2f%%)\n', ...
    theta_est(end, 2), A_true(1,2), 100 * abs(theta_est(end, 2) - A_true(1,2))/abs(A_true(1,2)));
fprintf('a21: %.4f (True: %.4f - Error: %.2f%%)\n', ...
    theta_est(end, 3), A_true(2,1), 100 * abs(theta_est(end, 3) - A_true(2,1))/abs(A_true(2,1)));
fprintf('a22: %.4f (True: %.4f - Error: %.2f%%)\n', ...
    theta_est(end, 4), A_true(2,2), 100 * abs(theta_est(end, 4) - A_true(2,2))/abs(A_true(2,2)));
fprintf('b1:  %.4f (True: %.4f - Error: %.2f%%)\n', ...
    theta_est(end, 5), B_true(1), 100 * abs(theta_est(end, 5) - B_true(1)));
fprintf('b2:  %.4f (True: %.4f - Error: %.2f%%)\n', ...
    theta_est(end, 6), B_true(2), 100 * abs(theta_est(end, 6) - B_true(2))/abs(B_true(2)));

fprintf('\nFinal estimation error norm: %.6f\n', sqrt(e_x(end,1)^2 + e_x(end,2)^2));
fprintf('Final parameter norm: %.4f (Bound M = %.1f)\n', theta_norms(end), M);
