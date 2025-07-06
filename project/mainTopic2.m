%% Modeling and Simulation of Dynamic Systems - Project - Topic 2
% Author: Aristeidis Daskalopoulos - AEM:10640

clear; clc; close all;
addpath('src\');
addpath('utils\');

%% 1. Define the True System Parameters
theta1_true = 1.5;  % Parameter θ1
theta2_true = 2.0;  % Parameter θ2

% True system function
f_true = @(x, u) -x^3 + theta1_true * tanh(x) + theta2_true / (1 + x^2) + u;

%% 2. Design Input Signal u(t)
% Rich input signal - sum of sinusoids
input_signal = @(t) sin(0.5*t) + 0.8*sin(1.2*t) + 0.6*sin(2.3*t) + 0.4*sin(3.7*t);

%% 3. Simulate True System
T_final = 100;  % Total simulation time
dt = 0.1;       % Time step
t  = 0:dt:T_final;
N  = length(t);

% Generate input data
u_data = zeros(N, 1);
for i = 1:N
    u_data(i) = input_signal(t(i));
end

% Simulate true system using ODE solver
x0_true = 0;  % Initial condition
[t_sim, x_true] = ode45(@(t, x) f_true(x, input_signal(t)), [0 T_final], x0_true);

% Interpolate to get uniform sampling
x_true_uniform = interp1(t_sim, x_true, t);

%% 4. Filter Design Parameters
lambda = 1.0;

%% 5. Model Structure Definitions

% Model 1: Simple Polynomial
model1_basis = @(x, u) [1; x; x^2; x^3; u];  % 5 parameters

% Model 2: Polynomial + Tanh
model2_basis = @(x, u) [1; x; x^3; tanh(x); u];  % 5 parameters

% Model 3: Gaussian RBFs
rbf_centers = [-2, -1, 0, 1, 2];  % RBF centers
rbf_sigma = 1.0;  % RBF width
model3_basis = @(x, u) [exp(-(x - rbf_centers(1))^2/(2*rbf_sigma^2)); 
                        exp(-(x - rbf_centers(2))^2/(2*rbf_sigma^2));
                        exp(-(x - rbf_centers(3))^2/(2*rbf_sigma^2));
                        exp(-(x - rbf_centers(4))^2/(2*rbf_sigma^2));
                        exp(-(x - rbf_centers(5))^2/(2*rbf_sigma^2));
                        u];  % 6 parameters

%% 6. Cross-Validation Setup
num_folds   =  4;
fold_size   =  floor(N / num_folds);
models      =  {model1_basis, model2_basis, model3_basis};
model_names =  {'Polynomial', 'Polynomial + Tanh', 'Gaussian RBFs'};
num_models  =  length(models);

% Initialize results storage
cv_mse    =  zeros(num_models, num_folds);
cv_params =  cell(num_models, num_folds);

%% 7. Cross-Validation Loop
fprintf('Starting Cross-Validation...\n');

for model_idx = 1:num_models
    fprintf('Evaluating %s model...\n', model_names{model_idx});
    
    for fold = 1:num_folds
        % Define test and training indices
        test_start   = (fold-1) * fold_size + 1;
        test_end     = min(fold * fold_size, N);
        test_indices = test_start:test_end;
        
        train_indices = setdiff(1:N, test_indices);
        
        % Training data
        t_train = t(train_indices);
        x_train = x_true_uniform(train_indices);
        u_train = u_data(train_indices);
        
        % Test data
        t_test = t(test_indices);
        x_test = x_true_uniform(test_indices);
        u_test = u_data(test_indices);
        
        % % Debug: Check data sizes
        % fprintf('  Fold %d: Train=%d, Test=%d samples\n', ...
        %     fold, length(train_indices), length(test_indices));
        
        % Parameter estimation for current model
        try
            [theta_est, x_hat_train] = ...
                estimateParameters(t_train, x_train, u_train, models{model_idx}, lambda);
            
            % Evaluate on test set
            x_hat_test = ...
                evaluateModel(t_test, x_test, u_test, theta_est, models{model_idx}, lambda);
            
            % Ensure both are column vectors for proper subtraction
            if size(x_test, 2) > size(x_test, 1)
                x_test = x_test';
            end
            if size(x_hat_test, 2) > size(x_hat_test, 1)
                x_hat_test = x_hat_test';
            end
            
            % % Debug: Check final sizes
            % fprintf('    x_test size: [%d, %d], x_hat_test size: [%d, %d]\n', ...
            %         size(x_test, 1), size(x_test, 2), size(x_hat_test, 1), size(x_hat_test, 2));
            
            % Compute MSE for this fold
            cv_mse(model_idx, fold)    = mean((x_test - x_hat_test).^2);
            cv_params{model_idx, fold} = theta_est;
            
        catch ME
            fprintf('    Error in fold %d: %s\n', fold, ME.message);
            cv_mse(model_idx, fold)    = inf;  % Set to infinity for failed folds :(
            cv_params{model_idx, fold} = [];
        end
    end
end

%% 8. Compute Final Metrics
mean_cv_mse = mean(cv_mse, 2);
std_cv_mse  = std(cv_mse, 0, 2);

% Train on full dataset for final evaluation
fprintf('\nTraining final models on full dataset...\n');
final_params = cell(num_models, 1);
final_x_hat  = zeros(N, num_models);
final_mse    = zeros(num_models, 1);
aic_scores   = zeros(num_models, 1);

for model_idx = 1:num_models
    [theta_final, x_hat_final] = ...
        estimateParameters(t, x_true_uniform, u_data, models{model_idx}, lambda);

    % Ensure both are column vectors for proper subtraction
    if size(x_hat_final, 2) > size(x_hat_final, 1)
        x_hat_final = x_hat_final';
    end
    if size(x_true_uniform, 2) > size(x_true_uniform, 1)
        x_true_uniform = x_true_uniform';
    end
    
    final_params{model_idx}   = theta_final;
    final_x_hat(:, model_idx) = x_hat_final;
    final_mse(model_idx)      = mean((x_true_uniform - x_hat_final).^2);
    
    % Compute AIC
    num_params            = length(theta_final);
    aic_scores(model_idx) = 2 * num_params + N * log(final_mse(model_idx));
end

%% 9. Display Results
fprintf('\n==== MODEL EVALUATION RESULTS ====\n');
fprintf('%-20s | %-12s | %-12s | %-12s | %-12s\n', ...
        'Model', 'CV MSE', 'CV Std', 'Final MSE', 'AIC');
fprintf('%-20s-|%-12s-|%-12s-|%-12s-|%-12s\n', ...
        repmat('-', 1, 20), repmat('-', 1, 12), repmat('-', 1, 12), ...
        repmat('-', 1, 12), repmat('-', 1, 12));

for i = 1:num_models
    fprintf('%-20s | %-12.6f | %-12.6f | %-12.6f | %-12.2f\n', ...
            model_names{i}, mean_cv_mse(i), std_cv_mse(i), final_mse(i), aic_scores(i));
end

% Find best model
[~, best_model_idx] = min(mean_cv_mse);
fprintf('\nBest model based on CV MSE: %s\n', model_names{best_model_idx});

%% 10. Plotting Results
% Plot 1: True vs Estimated trajectories
figure('Name', 'System Trajectories Comparison', 'NumberTitle', 'off');
subplot(2, 2, 1);
plot(t, x_true_uniform, 'b-', 'LineWidth', 2);
hold on;
plot(t, final_x_hat(:, 1), 'r--', 'LineWidth', 1.5);
plot(t, final_x_hat(:, 2), 'g--', 'LineWidth', 1.5);
plot(t, final_x_hat(:, 3), 'm--', 'LineWidth', 1.5);
title('System Output Comparison');
xlabel('Time [s]');
ylabel('x(t)');
legend('True', model_names{:}, 'Location', 'best');
grid on;

% Plot 2: Input signal
subplot(2, 2, 2);
plot(t, u_data, 'k-', 'LineWidth', 1.5);
title('Input Signal u(t)');
xlabel('Time [s]');
ylabel('u(t)');
grid on;

% Plot 3: Estimation errors
subplot(2, 2, 3);
for i = 1:num_models
    plot(t, x_true_uniform - final_x_hat(:, i), 'LineWidth', 1.5);
    hold on;
end
title('Estimation Errors');
xlabel('Time [s]');
ylabel('Error = x_{true} - x_{hat}');
legend(model_names{:}, 'Location', 'best');
grid on;

% Plot 4: Cross-validation results
subplot(2, 2, 4);
bar(mean_cv_mse);
hold on;
errorbar(1:num_models, mean_cv_mse, std_cv_mse, 'k.', 'LineWidth', 1.5);
title('Cross-Validation MSE');
xlabel('Model');
ylabel('MSE');
set(gca, 'XTickLabel', model_names);
grid on;

% Plot 5: Model comparison metrics
figure('Name', 'Model Comparison Metrics', 'NumberTitle', 'off');
subplot(1, 2, 1);
bar([mean_cv_mse, final_mse]);
title('MSE Comparison');
xlabel('Model');
ylabel('MSE');
legend('CV MSE', 'Final MSE', 'Location', 'best');
set(gca, 'XTickLabel', model_names);
grid on;

subplot(1, 2, 2);
bar(aic_scores);
title('AIC Scores (Lower -> Better)');
xlabel('Model');
ylabel('AIC');
set(gca, 'XTickLabel', model_names);
grid on;

%% 11. Parameter Analysis for Best Model
fprintf('\n==== BEST MODEL PARAMETER ANALYSIS ====\n');
fprintf('Best Model: %s\n', model_names{best_model_idx});
fprintf('Final Parameters:\n');
best_params = final_params{best_model_idx};
for i = 1:length(best_params)
    fprintf('θ_%d = %.4f\n', i, best_params(i));
end

% Additional analysis plots
figure('Name', 'Best Model Analysis', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t, x_true_uniform, 'b-', 'LineWidth', 2);
hold on;
plot(t, final_x_hat(:, best_model_idx), 'r--', 'LineWidth', 1.5);
title(sprintf('Best Model (%s) vs True System', model_names{best_model_idx}));
xlabel('Time [s]');
ylabel('x(t)');
legend('True System', 'Best Model', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, x_true_uniform - final_x_hat(:, best_model_idx), 'g-', 'LineWidth', 1.5);
title('Best Model Estimation Error');
xlabel('Time [s]');
ylabel('Error');
grid on;

% Calculate final error statistics
final_error =  x_true_uniform - final_x_hat(:, best_model_idx);
rmse        =  sqrt(mean(final_error.^2));
mae         =  mean(abs(final_error));

fprintf('\nFinal Error Statistics for Best Model:\n');
fprintf('RMSE: %.6f\n', rmse);
fprintf('MAE:  %.6f\n', mae);
