%% Nonlinear System Model Selection and Evaluation
% This code implements model selection and evaluation for approximating
% an unknown nonlinear dynamic system using different basis functions

clear all; close all; clc;

%% 1. TRUE SYSTEM PARAMETERS AND SETUP
% Choose parameters within the specified range [0.5, 2]
theta1_true = 1.2;  % True parameter 1
theta2_true = 0.8;  % True parameter 2
theta_true = [theta1_true; theta2_true];

% Time vector
t_span = [0 10];  % 10 seconds simulation
dt = 0.01;
t = 0:dt:10;

% Input signal u(t) - choosing a rich excitation signal
u_func = @(t) 0.5*sin(2*pi*0.1*t) + 0.3*sin(2*pi*0.3*t) + 0.2*sin(2*pi*0.7*t);

% True system dynamics
f_true = @(x, u, theta) -x^3 + theta(1)*tanh(x) + theta(2)/(1+x^2) + u;

% Generate true system data
x0 = 0.1;  % Initial condition
[t_true, x_true] = ode45(@(t,x) f_true(x, u_func(t), theta_true), t_span, x0);

% Add measurement noise to make it more realistic
noise_level = 0.01;
x_measured = x_true + noise_level * randn(size(x_true));

% Calculate input values at time points
u_values = u_func(t_true);

%% 2. DEFINE CANDIDATE MODEL STRUCTURES
% Each model structure uses different basis functions to approximate f(x,u,θ)

% Model 1: Polynomial basis
% f_model1 = θ1*x + θ2*x^2 + θ3*x^3 + u
phi1 = @(x, u) [x, x.^2, x.^3, ones(size(x))];  % Last term for u coefficient

% Model 2: Polynomial + Trigonometric basis
% f_model2 = θ1*x + θ2*x^2 + θ3*sin(x) + θ4*cos(x) + u
phi2 = @(x, u) [x, x.^2, sin(x), cos(x), ones(size(x))];

% Model 3: Polynomial + Hyperbolic basis (closer to true system)
% f_model3 = θ1*x + θ2*x^3 + θ3*tanh(x) + θ4/(1+x^2) + u
phi3 = @(x, u) [x, x.^3, tanh(x), 1./(1+x.^2), ones(size(x))];

% Model 4: Radial Basis Functions (Gaussian)
% f_model4 = θ1*exp(-(x-c1)^2) + θ2*exp(-(x-c2)^2) + θ3*exp(-(x-c3)^2) + u
c1 = -1; c2 = 0; c3 = 1;  % Centers for RBF
phi4 = @(x, u) [exp(-(x-c1).^2), exp(-(x-c2).^2), exp(-(x-c3).^2), ones(size(x))];

% Model 5: Extended polynomial basis
% f_model5 = θ1*x + θ2*x^2 + θ3*x^3 + θ4*x^4 + θ5*x^5 + u
phi5 = @(x, u) [x, x.^2, x.^3, x.^4, x.^5, ones(size(x))];

% Store all models
models = {phi1, phi2, phi3, phi4, phi5};
model_names = {'Polynomial', 'Poly+Trig', 'Poly+Hyperbolic', 'RBF', 'Extended Poly'};
num_models = length(models);

%% 3. PARAMETER ESTIMATION METHODS
% We'll use both batch least squares and recursive least squares

%% 3.1 Batch Least Squares Estimation
fprintf('=== BATCH LEAST SQUARES ESTIMATION ===\n');

% Calculate derivatives using finite differences
x_dot = gradient(x_measured, dt);

model_params = cell(num_models, 1);
model_errors = zeros(num_models, 1);
model_predictions = cell(num_models, 1);

for i = 1:num_models
    % Build regression matrix
    Phi = models{i}(x_measured, u_values);
    
    % Least squares estimation: θ = (Φ^T Φ)^(-1) Φ^T x_dot
    theta_est = pinv(Phi) * x_dot;
    model_params{i} = theta_est;
    
    % Calculate model prediction
    x_dot_pred = Phi * theta_est;
    model_predictions{i} = x_dot_pred;
    
    % Calculate mean squared error
    mse = mean((x_dot - x_dot_pred).^2);
    model_errors(i) = mse;
    
    fprintf('Model %d (%s): MSE = %.6f\n', i, model_names{i}, mse);
    fprintf('  Parameters: [%.4f', theta_est(1));
    for j = 2:length(theta_est)
        fprintf(', %.4f', theta_est(j));
    end
    fprintf(']\n');
end

%% 3.2 Recursive Least Squares (RLS) Implementation
fprintf('\n=== RECURSIVE LEAST SQUARES ESTIMATION ===\n');

% Initialize RLS parameters
forgetting_factor = 0.99;  % Forgetting factor for RLS
rls_results = cell(num_models, 1);

for i = 1:num_models
    % Get dimension of parameter vector
    phi_sample = models{i}(x_measured(1), u_values(1));
    n_params = length(phi_sample);
    
    % Initialize RLS
    theta_rls = zeros(n_params, length(t_true));
    P = 1000 * eye(n_params);  % Initial covariance matrix
    theta_curr = zeros(n_params, 1);  % Initial parameter estimate
    
    % RLS recursion
    for k = 2:length(t_true)
        % Get regression vector
        phi_k = models{i}(x_measured(k), u_values(k))';
        
        % Update covariance matrix
        P = (P - (P * (phi_k * phi_k') * P) / (forgetting_factor + phi_k' * P * phi_k)) / forgetting_factor;
        
        % Update parameter estimate
        K = P * phi_k / (forgetting_factor + phi_k' * P * phi_k);  % Kalman gain
        theta_curr = theta_curr + K * (x_dot(k) - phi_k' * theta_curr);
        
        theta_rls(:, k) = theta_curr;
    end
    
    rls_results{i} = theta_rls;
    
    % Final RLS error
    Phi = models{i}(x_measured, u_values);
    x_dot_pred_rls = Phi * theta_curr;
    mse_rls = mean((x_dot - x_dot_pred_rls).^2);
    
    fprintf('Model %d (%s) RLS: Final MSE = %.6f\n', i, model_names{i}, mse_rls);
end

%% 4. MODEL VALIDATION AND SIMULATION
fprintf('\n=== MODEL VALIDATION THROUGH SIMULATION ===\n');

% Generate validation data with different input
u_val_func = @(t) 0.8*sin(2*pi*0.2*t) + 0.4*cos(2*pi*0.5*t);
[t_val, x_val_true] = ode45(@(t,x) f_true(x, u_val_func(t), theta_true), t_span, x0);

simulation_errors = zeros(num_models, 1);

for i = 1:num_models
    % Create model ODE function
    theta_model = model_params{i};
    
    % Model dynamics: x_dot = phi(x,u) * theta
    model_ode = @(t, x) models{i}(x, u_val_func(t)) * theta_model;
    
    % Simulate model
    [t_sim, x_sim] = ode45(model_ode, t_span, x0);
    
    % Interpolate to compare at same time points
    x_sim_interp = interp1(t_sim, x_sim, t_val, 'linear', 'extrap');
    
    % Calculate simulation error
    sim_error = sqrt(mean((x_val_true - x_sim_interp).^2));
    simulation_errors(i) = sim_error;
    
    fprintf('Model %d (%s): Simulation RMSE = %.6f\n', i, model_names{i}, sim_error);
end

%% 5. INFORMATION CRITERIA FOR MODEL SELECTION
fprintf('\n=== INFORMATION CRITERIA ===\n');

n_data = length(x_measured);
aic_values = zeros(num_models, 1);
bic_values = zeros(num_models, 1);

for i = 1:num_models
    n_params = length(model_params{i});
    mse = model_errors(i);
    
    % AIC = n*log(MSE) + 2*k
    aic_values(i) = n_data * log(mse) + 2 * n_params;
    
    % BIC = n*log(MSE) + k*log(n)
    bic_values(i) = n_data * log(mse) + n_params * log(n_data);
    
    fprintf('Model %d (%s): AIC = %.2f, BIC = %.2f\n', i, model_names{i}, aic_values(i), bic_values(i));
end

%% 6. CROSS-VALIDATION
fprintf('\n=== CROSS-VALIDATION ===\n');

k_fold = 5;  % 5-fold cross-validation
cv_errors = zeros(num_models, k_fold);

% Divide data into k folds
fold_size = floor(length(x_measured) / k_fold);

for fold = 1:k_fold
    % Define test indices
    test_start = (fold-1) * fold_size + 1;
    test_end = min(fold * fold_size, length(x_measured));
    test_idx = test_start:test_end;
    
    % Training indices (everything except test)
    train_idx = setdiff(1:length(x_measured), test_idx);
    
    % Training data
    x_train = x_measured(train_idx);
    u_train = u_values(train_idx);
    x_dot_train = x_dot(train_idx);
    
    % Test data
    x_test = x_measured(test_idx);
    u_test = u_values(test_idx);
    x_dot_test = x_dot(test_idx);
    
    % Train and test each model
    for i = 1:num_models
        % Train on training data
        Phi_train = models{i}(x_train, u_train);
        theta_fold = pinv(Phi_train) * x_dot_train;
        
        % Test on test data
        Phi_test = models{i}(x_test, u_test);
        x_dot_pred_test = Phi_test * theta_fold;
        
        % Calculate test error
        cv_errors(i, fold) = mean((x_dot_test - x_dot_pred_test).^2);
    end
end

cv_mean_errors = mean(cv_errors, 2);
cv_std_errors = std(cv_errors, 0, 2);

for i = 1:num_models
    fprintf('Model %d (%s): CV Error = %.6f ± %.6f\n', i, model_names{i}, cv_mean_errors(i), cv_std_errors(i));
end

%% 7. VISUALIZATION
% Create comprehensive plots

figure('Position', [100, 100, 1200, 800]);

% Plot 1: True system response
subplot(2, 3, 1);
plot(t_true, x_true, 'b-', 'LineWidth', 2);
hold on;
plot(t_true, x_measured, 'r.', 'MarkerSize', 4);
xlabel('Time (s)');
ylabel('x(t)');
title('True System Response');
legend('True x(t)', 'Measured x(t)', 'Location', 'best');
grid on;

% Plot 2: Model comparison - derivative prediction
subplot(2, 3, 2);
plot(t_true, x_dot, 'b-', 'LineWidth', 2);
hold on;
colors = {'r--', 'g--', 'm--', 'c--', 'k--'};
for i = 1:num_models
    plot(t_true, model_predictions{i}, colors{i}, 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('dx/dt');
title('Derivative Prediction Comparison');
legend(['True', model_names], 'Location', 'best');
grid on;

% Plot 3: Model errors comparison
subplot(2, 3, 3);
bar(1:num_models, model_errors);
set(gca, 'XTickLabel', model_names);
ylabel('Mean Squared Error');
title('Model Fitting Errors');
xtickangle(45);
grid on;

% Plot 4: Information criteria
subplot(2, 3, 4);
bar(1:num_models, [aic_values, bic_values]);
set(gca, 'XTickLabel', model_names);
ylabel('Information Criterion Value');
title('AIC and BIC Comparison');
legend('AIC', 'BIC', 'Location', 'best');
xtickangle(45);
grid on;

% Plot 5: Cross-validation results
subplot(2, 3, 5);
errorbar(1:num_models, cv_mean_errors, cv_std_errors, 'bo-', 'LineWidth', 2);
set(gca, 'XTickLabel', model_names);
ylabel('CV Error');
title('Cross-Validation Results');
xtickangle(45);
grid on;

% Plot 6: Simulation validation
subplot(2, 3, 6);
bar(1:num_models, simulation_errors);
set(gca, 'XTickLabel', model_names);
ylabel('Simulation RMSE');
title('Simulation Validation Errors');
xtickangle(45);
grid on;

%% 8. RLS PARAMETER EVOLUTION
figure('Position', [100, 100, 1200, 600]);

% Plot parameter evolution for the best model (Model 3 - closest to true structure)
best_model_idx = 3;  % Poly+Hyperbolic should be closest to true system

subplot(1, 2, 1);
plot(t_true, rls_results{best_model_idx}');
xlabel('Time (s)');
ylabel('Parameter Value');
title(sprintf('RLS Parameter Evolution - %s', model_names{best_model_idx}));
legend('θ₁', 'θ₂', 'θ₃', 'θ₄', 'θ₅', 'Location', 'best');
grid on;

% Plot parameter convergence for all models
subplot(1, 2, 2);
final_params = zeros(num_models, 1);
for i = 1:num_models
    final_params(i) = norm(rls_results{i}(:, end) - rls_results{i}(:, end-100));
end
bar(1:num_models, final_params);
set(gca, 'XTickLabel', model_names);
ylabel('Parameter Change (Last 100 samples)');
title('Parameter Convergence');
xtickangle(45);
grid on;

%% 9. SUMMARY AND RECOMMENDATIONS
fprintf('\n=== SUMMARY AND RECOMMENDATIONS ===\n');
fprintf('Model Selection Results:\n');
fprintf('1. Best fitting model (MSE): %s (MSE = %.6f)\n', model_names{find(model_errors == min(model_errors))}, min(model_errors));
fprintf('2. Best model (AIC): %s (AIC = %.2f)\n', model_names{find(aic_values == min(aic_values))}, min(aic_values));
fprintf('3. Best model (BIC): %s (BIC = %.2f)\n', model_names{find(bic_values == min(bic_values))}, min(bic_values));
fprintf('4. Best model (CV): %s (CV Error = %.6f)\n', model_names{find(cv_mean_errors == min(cv_mean_errors))}, min(cv_mean_errors));
fprintf('5. Best model (Simulation): %s (RMSE = %.6f)\n', model_names{find(simulation_errors == min(simulation_errors))}, min(simulation_errors));

% Overall ranking
overall_ranks = zeros(num_models, 1);
[~, rank_mse] = sort(model_errors);
[~, rank_aic] = sort(aic_values);
[~, rank_bic] = sort(bic_values);
[~, rank_cv] = sort(cv_mean_errors);
[~, rank_sim] = sort(simulation_errors);

for i = 1:num_models
    overall_ranks(i) = find(rank_mse == i) + find(rank_aic == i) + find(rank_bic == i) + find(rank_cv == i) + find(rank_sim == i);
end

[~, best_overall] = min(overall_ranks);
fprintf('\nOverall best model: %s\n', model_names{best_overall});
fprintf('True system parameters: θ₁ = %.3f, θ₂ = %.3f\n', theta1_true, theta2_true);

if best_overall == 3  % If Poly+Hyperbolic is selected
    fprintf('Estimated parameters: θ₁ = %.3f, θ₂ = %.3f\n', model_params{3}(3), model_params{3}(4));
    fprintf('Parameter estimation error: Δθ₁ = %.3f, Δθ₂ = %.3f\n', ...
        abs(model_params{3}(3) - theta1_true), abs(model_params{3}(4) - theta2_true));
end

fprintf('\nThis analysis demonstrates the importance of:\n');
fprintf('1. Using multiple evaluation criteria\n');
fprintf('2. Validating models through simulation\n');
fprintf('3. Considering model complexity vs. accuracy trade-offs\n');
fprintf('4. Using proper cross-validation techniques\n');
