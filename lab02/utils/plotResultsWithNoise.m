function plotResultsWithNoise(t, x, x_hat, x_measured, e_x, m_est, b_est, k_est, m_true, b_true, k_true, T_final, titleStr)
    % Figure 1: States and Error
    figure('Name', [titleStr ' - States']);
    
    % Position plot with noisy measurement
    subplot(2,1,1);
    plot(t, x_measured, 'g:', t, x, 'b-', t, x_hat, 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Position x(t) [m]');
    title('Position: x(t), x_{measured}(t) and x̂(t)');
    legend('Measured x(t)', 'True x(t)', 'Estimated x̂(t)');
    
    % Error plot
    subplot(2,1,2);
    plot(t, e_x, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Error e_x(t) [m]');
    title('Position Error: e_x(t) = x_{measured}(t) - x̂(t)');
    
    % Parameter Estimates
    figure('Name', [titleStr ' - Parameters']);
    plot(t, m_est, 'b-', ...
         t, b_est, 'r--', ...
         t, k_est, 'g:', ...
         [0 T_final], [m_true m_true], 'b--', ...
         [0 T_final], [b_true b_true], 'r:', ...
         [0 T_final], [k_true k_true], 'g-.', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Parameter Values');
    title('Parameter Estimates');
    legend('m̂(t)', 'b̂(t)', 'k̂(t)', 'm_{true}', 'b_{true}', 'k_{true}');
end