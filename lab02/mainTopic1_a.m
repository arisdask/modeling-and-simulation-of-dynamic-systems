%% Modeling and Simulation of Dynamic Systems - Assignment 2
% Author: Aristeidis Daskalopoulos - AEM:10640

% Mass-Spring-Damper System - Parameter Estimation using Gradient Method
% This script implements the gradient method to estimate the parameters 
% of a mass-spring-damper system in real-time (on-line).
%
% System:     m*d²x/dt² + b*dx/dt + k*x = u(t)
% True values: m = 1.315, b = 0.225, k = 0.725

clear; clc; close all;

addpath('src');
addpath('utils');


%% System Parameters
m_true = 1.315;    % True mass value
b_true = 0.225;    % True damping coefficient
k_true = 0.725;    % True spring stiffness

% Choose filter version (1: first-order filter, 2: second-order filter)
filter_version = 2;

% Filter parameters
lambda1 = 0.3;     % First filter parameter (must be positive for stability)
lambda2 = 0.02;    % Second filter parameter (for version 2 only)

gamma = 1;     % Learning rate

% Simulation parameters
T_final = 20;     % Simulation time (seconds)
tspan = [0 T_final]; % Time span

% Set up initial conditions based on filter version
if filter_version == 1
    % Version 1: First-order filter
    x0 = [0; 0];              % Initial state [position; velocity]
    phi0 = [0; 0; 0];         % Initial filtered signals
    theta0 = [0; 0.5; 0.5]; % Initial parameter estimates: [a1-λ1; a2; b0]
    initial_conditions = [x0; phi0; theta0];
else
    % Version 2: Second-order filter
    x0 = [0; 0];              % Initial state [position; velocity]
    phi10 = [0; 0];           % Initial phi1 and its derivative
    phi20 = [0; 0];           % Initial phi2 and its derivative
    phi30 = [0; 0];           % Initial phi3 and its derivative
    theta0 = [-0.2; 0.48; 1]; % Initial parameter estimates: [a1-λ1; a2-λ2; b0]
    initial_conditions = [x0; phi10; phi20; phi30; theta0];
end


%% Case (i): Constant Input u(t) = 2.5
input_type = 1;  % Stands for constant input case

% Run simulation with appropriate dynamics function
if filter_version == 1
    [t1, y1] = ode45( ...
        @(t, y) systemDynamics1(t, y, m_true, b_true, k_true, lambda1, gamma, input_type), ...
        tspan, ...
        initial_conditions ...
        );
    
    % Extract and process results for version 1
    [x_case1, xdot_case1, phi_case1, theta_case1, ...
        m_est_case1, b_est_case1, k_est_case1, ...
        x_hat_case1, xdot_hat_case1, e_x_case1, e_xdot_case1] = extractResults1(t1, y1, lambda1);
else
    [t1, y1] = ode45( ...
        @(t, y) systemDynamics2(t, y, m_true, b_true, k_true, lambda1, lambda2, gamma, input_type), ...
        tspan, ...
        initial_conditions ...
        );
    
    % Extract and process results for version 2
    [x_case1, xdot_case1, phi_case1, theta_case1, ...
        m_est_case1, b_est_case1, k_est_case1, ...
        x_hat_case1, xdot_hat_case1, e_x_case1, e_xdot_case1] = extractResults2(t1, y1, lambda1, lambda2);
end

% Display final parameter estimates for case (i)
if filter_version == 1
    title_str = sprintf('Constant Input u(t) = 2.5 | gamma = %.3f | lambda = %.2f', gamma, lambda1);
else
    title_str = sprintf('Constant Input u(t) = 2.5 | gamma = %.3f | lambda1 = %.2f, lambda2 = %.2f', ...
        gamma, lambda1, lambda2);
end

displayResults( ...
    m_true, b_true, k_true, ...
    mean(m_est_case1(end-5:end)), mean(b_est_case1(end-5:end)), mean(k_est_case1(end-5:end)), ...
    title_str ...
    );

% Plot results for case (i)
plotResults( ...
    t1, x_case1, xdot_case1, x_hat_case1, xdot_hat_case1, e_x_case1, e_xdot_case1, ...
    m_est_case1, b_est_case1, k_est_case1, m_true, b_true, k_true, ...
    title_str ...
    );


%% Case (ii): Sinusoidal Input u(t) = 2.5*sin(t)
input_type = 2;  % Stands for sinusoidal input case

% Run simulation with appropriate dynamics function
if filter_version == 1
    [t2, y2] = ode45( ...
        @(t, y) systemDynamics1(t, y, m_true, b_true, k_true, lambda1, gamma, input_type), ...
        tspan, ...
        initial_conditions ...
        );
    
    % Extract and process results for version 1
    [x_case2, xdot_case2, phi_case2, theta_case2, ...
        m_est_case2, b_est_case2, k_est_case2, ...
        x_hat_case2, xdot_hat_case2, e_x_case2, e_xdot_case2] = extractResults1(t2, y2, lambda1);
else
    [t2, y2] = ode45( ...
        @(t, y) systemDynamics2(t, y, m_true, b_true, k_true, lambda1, lambda2, gamma, input_type), ...
        tspan, ...
        initial_conditions ...
        );
    
    % Extract and process results for version 2
    [x_case2, xdot_case2, phi_case2, theta_case2, ...
        m_est_case2, b_est_case2, k_est_case2, ...
        x_hat_case2, xdot_hat_case2, e_x_case2, e_xdot_case2] = extractResults2(t2, y2, lambda1, lambda2);
end

% Display final parameter estimates for case (ii)
if filter_version == 1
    title_str = sprintf('Sinusoidal Input u(t) = 2.5*sin(t) | gamma = %.3f | lambda = %.2f', gamma, lambda1);
else
    title_str = sprintf('Sinusoidal Input u(t) = 2.5*sin(t) | gamma = %.3f | lambda1 = %.2f, lambda2 = %.2f', gamma, lambda1, lambda2);
end

displayResults( ...
    m_true, b_true, k_true, ...
    mean(m_est_case2(end-5:end)), mean(b_est_case2(end-5:end)), mean(k_est_case2(end-5:end)), ...
    title_str ...
    );

% Plot results for case (ii)
plotResults( ...
    t2, x_case2, xdot_case2, x_hat_case2, xdot_hat_case2, e_x_case2, e_xdot_case2, ...
    m_est_case2, b_est_case2, k_est_case2, m_true, b_true, k_true, ...
    title_str ...
    );
