%% Vector Norm Calculation Script
clear; clc;
syms t

omega_bar = 0.3;

omega1 = 0.8 * omega_bar * sin(0.7 * t);
omega2 = 0.6 * omega_bar * cos(1.2 * t + pi/4);

omega_vector = [omega1; omega2];

fprintf('Original vector omega(t):\n');
disp(omega_vector);

norm_squared            =  omega1^2 + omega2^2;
norm_squared_expanded   =  expand(norm_squared);
norm_squared_simplified =  simplify(norm_squared_expanded);

fprintf('\nSimplified norm squared:\n');
disp(norm_squared_simplified);

norm_omega             =  sqrt(norm_squared_simplified);
norm_omega_simplified  =  simplify(norm_omega);

fprintf('\nSimplified norm ||omega(t)||:\n');
disp(norm_omega_simplified);


t_values    = linspace(0, 100, 1000);
norm_values = double(subs(norm_omega_simplified, t, t_values));

figure;
plot(t_values, norm_values, 'b-', 'LineWidth', 1.5);
xlabel('t');
ylabel('||ω(t)||');
title('Norm of Vector ω(t) vs Time');
grid on;
xlim([0, 100]);
ylim([0, max(norm_values) * 1.1]);
