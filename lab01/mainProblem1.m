% Modeling and Simulation of Dynamic Systems - Assignment 1
% Estimation of Unknown Parameters - Least Squares Method
% Author: Aristeidis Daskalopoulos - AEM:10640

%% Initialize Problem
clear; clc; close all;

addpath('src');
addpath('utils');

% System parameters
m = 0.75; L = 1.25; c = 0.15; g = 9.81;
A0 = 4; omega = 2;

% Input Torque u(t) function
inputTorque = @(t) A0 * sin(omega * t);

% Pendulum ODE function (for state x = [x1 x2]^T)
pendulumODE = @(t, x, m, L, c, g, inputFunc) [
    x(2);                                           % q' = dq/dt = x1' = x2
    (inputFunc(t) - c*x(2) - m*g*L*x(1)) / (m*L^2)  % q" = x2'
];

%% Problem 1 - Simulation
T_final = 20;    % total simulation time [sec]
dt      = 1e-4;  % integration step
time    = 0:dt:T_final;

% Initial conditions
x0 = [0; 0];  % [q(0) q'(0)] = [0 0]

% Solve the ODE
[t_sol, x_sol] = solveODE(pendulumODE, time, x0, m, L, c, g, inputTorque);

% Plot the results
plotPendulumStates(t_sol, x_sol);
plotInputTorque(time, inputTorque);
