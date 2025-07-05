%% Constrained Parameter Estimation for Linear Systems
% Implementation of mixed structure parameter estimation with constraints
% Author: Based on mathematical analysis provided

clear; clc; close all;
addpath('src\');
addpath('utils\');

%% True System Parameters
A_true = [-2.15,  0.25;
          -0.75, -2.0];

B_true = [0; 1.5];

%% Estimator Parameters
% Mixed structure parameter (for matrix C)
c_coeff = 10;  % Coefficient for matrix C (10)
C = c_coeff * eye(2);  % C = c_coeff * I

% Learning rates (diagonal elements of Gamma matrix)
gamma = ones(6,1) .* [5; 2; 5; 2; 2; 2;];
Gamma = diag(gamma);

%% Constraint Parameters
% Constraints: -3 ≤ a11 ≤ -1 and b2 ≥ 1
% g1 = -a11 - 3 ≤ 0  (a11 ≥ -3)
% g2 = a11 + 1 ≤ 0   (a11 ≤ -1)  
% g3 = -b2 + 1 ≤ 0   (b2 ≥ 1)

% Constraint gradients
nabla_g1 = [-1; 0; 0; 0; 0; 0];  % ∇g1
nabla_g2 = [1; 0; 0; 0; 0; 0];   % ∇g2
nabla_g3 = [0; 0; 0; 0; 0; -1];  % ∇g3

%% Simulation Parameters
T_final = 30;  % Simulation time (seconds)
tspan = [0 T_final];

%% Initial Conditions
% True system initial state
x0 = [1; 0.5];

% Initial parameter estimates (must be feasible)
% theta = [a11, a12, a21, a22, b1, b2]'
theta0 = [-2.5; 0.1; -0.5; -1.8; 1; 1.2];  % Satisfies constraints

% Initial estimated state (same as true state for mixed structure)
x_hat0 = x0;

% Combined initial conditions: [x; x_hat; theta]
initial_conditions = [x0; x_hat0; theta0];

%% Run Simulation
[t, x_sim] = ode45(@(t, x) ...
    constrainedSystemDynamics( ...
    t, x, A_true, B_true, C, Gamma, nabla_g1, nabla_g2, nabla_g3), ...
    tspan, initial_conditions);

%% Extract Results
n_states = 2;
n_params = 6;

% Extract state trajectories
x_true    = x_sim(:, 1:n_states);
x_hat     = x_sim(:, (n_states+1):(2*n_states));
theta_est = x_sim(:, (2*n_states+1):(2*n_states+n_params));

% Compute estimation error
e_x = x_true - x_hat;

% Extract individual parameter estimates
a11_est = theta_est(:, 1);
a12_est = theta_est(:, 2);
a21_est = theta_est(:, 3);
a22_est = theta_est(:, 4);
b1_est  = theta_est(:, 5);
b2_est  = theta_est(:, 6);

%% Plot Results
% State trajectories and estimation errors
figure;
subplot(1,2,1);
plot(t, x_true(:,1), 'b-', t, x_hat(:,1), 'r--', 'LineWidth', 1.5);
title('State x_1(t)');
xlabel('Time [s]');
ylabel('x_1');
legend('True', 'Estimated', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(t, x_true(:,2), 'b-', t, x_hat(:,2), 'r--', 'LineWidth', 1.5);
title('State x_2(t)');
xlabel('Time [s]');
ylabel('x_2');
legend('True', 'Estimated', 'Location', 'best');
grid on;

figure;
subplot(1,3,1);
plot(t, e_x(:,1), 'g-', 'LineWidth', 1.5);
title('Estimation Error e_1(t)');
xlabel('Time [s]');
ylabel('$e_1 = x_1 - \hat{x}_1$', 'Interpreter', 'latex');
grid on;

subplot(1,3,2);
plot(t, e_x(:,2), 'g-', 'LineWidth', 1.5);
title('Estimation Error e_2(t)');
xlabel('Time [s]');
ylabel('$e_2 = x_2 - \hat{x}_2$', 'Interpreter', 'latex');
grid on;

subplot(1,3,3);
plot(t, sqrt(e_x(:,1).^2 + e_x(:,2).^2), 'm-', 'LineWidth', 1.5);
title('Estimation Error Norm ||e||');
xlabel('Time [s]');
ylabel('||e||');
grid on;

% Parameter estimates for matrix A
figure;
subplot(2,2,1);
plot(t, a11_est, 'r-', 'LineWidth', 1.5);
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
plot(t, a12_est, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(1,2)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter a_{12}');
xlabel('Time [s]');
ylabel('a_{12}');
legend('Estimated', 'True', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(t, a21_est, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(2,1)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter a_{21}');
xlabel('Time [s]');
ylabel('a_{21}');
legend('Estimated', 'True', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(t, a22_est, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, A_true(2,2)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter a_{22}');
xlabel('Time [s]');
ylabel('a_{22}');
legend('Estimated', 'True', 'Location', 'best');
grid on;

% Parameter estimates for matrix B
figure;
subplot(1,2,1);
plot(t, b1_est, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, B_true(1)*ones(size(t)), 'b--', 'LineWidth', 1.5);
title('Parameter b_1');
xlabel('Time [s]');
ylabel('b_1');
legend('Estimated', 'True', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(t, b2_est, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, B_true(2)*ones(size(t)), 'b--', 'LineWidth', 1.5);
plot(t, ones(size(t)), 'k:', 'LineWidth', 1);
title('Parameter b_2');
xlabel('Time [s]');
ylabel('b_2');
legend('Estimated', 'True', 'Constraint b_2 \geq 1', 'Location', 'best');
grid on;

%% Display final estimates
fprintf('Final Parameter Estimates:\n');
fprintf('a11: %.4f (True: %.4f)\n', a11_est(end), A_true(1,1));
fprintf('a12: %.4f (True: %.4f)\n', a12_est(end), A_true(1,2));
fprintf('a21: %.4f (True: %.4f)\n', a21_est(end), A_true(2,1));
fprintf('a22: %.4f (True: %.4f)\n', a22_est(end), A_true(2,2));
fprintf('b1:  %.4f (True: %.4f)\n', b1_est(end), B_true(1));
fprintf('b2:  %.4f (True: %.4f)\n', b2_est(end), B_true(2));

%% Check constraint satisfaction
% fprintf('\nConstraint Satisfaction:\n');
% fprintf('a11 ∈ [-3, -1]: %.4f ∈ [%.1f, %.1f] - %s\n', a11_est(end), -3, -1, ...
%     (a11_est(end) >= -3 && a11_est(end) <= -1) ? 'OK' : 'VIOLATED');
% fprintf('b2 ≥ 1: %.4f ≥ 1 - %s\n', b2_est(end), (b2_est(end) >= 1) ? 'OK' : 'VIOLATED');
