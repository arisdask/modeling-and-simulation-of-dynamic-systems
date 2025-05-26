%% Function to plot results
function plotResults(t, x, xdot, x_hat, xdot_hat, e_x, e_xdot, m_est, b_est, k_est, m_true, b_true, k_true, title_text)
    % Create figures with subplots
    
    % Figure 1: States and their errors
    figure;
    
    % Position state and estimate
    subplot(2, 2, 1);
    plot(t, x, 'b-', t, x_hat, 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Position x(t) [m]');
    title('Position: x(t) and x̂(t)');
    legend('True x(t)', 'Estimated x_{hat}(t)');
    
    % Position error
    subplot(2, 2, 2);
    plot(t, e_x, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Error e_x(t) [m]');
    title('Position Error: e_x(t) = x(t) - x_{hat}(t)');
    
    % Velocity state and estimate
    subplot(2, 2, 3);
    plot(t, xdot, 'b-', t, xdot_hat, 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel("Velocity x'(t) [m/s]");
    title("Velocity: x'(t) and x'_{hat}(t)");
    legend("True x'(t)", "Estimated x'_{hat}(t)");
    
    % Velocity error
    subplot(2, 2, 4);
    plot(t, e_xdot, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Error [m/s]');
    title("Velocity Error: e_x'(t) = x'(t) - x'_{hat}(t)");
    
    % Add input type as suptitle
    sgtitle(['System States and Errors - ' title_text], 'FontSize', 14);
    
    % Figure 2: Parameter estimates
    figure;
    
    % Mass estimate
    subplot(1, 3, 1);
    plot(t, m_est, 'b-', [t(1), t(end)], [m_true, m_true], 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Mass (kg)');
    title('Mass Estimate: m̂(t)');
    legend('Estimate', 'True Value');
    
    % Damping coefficient estimate
    subplot(1, 3, 2);
    plot(t, b_est, 'b-', [t(1), t(end)], [b_true, b_true], 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Damping Coef.');
    title('Damping Coefficient Estimate: b̂(t)');
    legend('Estimate', 'True Value');
    
    % Spring constant estimate
    subplot(1, 3, 3);
    plot(t, k_est, 'b-', [t(1), t(end)], [k_true, k_true], 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Spring Const.');
    title('Spring Constant Estimate: k̂(t)');
    legend('Estimate', 'True Value');
    
    % Add input type as suptitle
    sgtitle(['Parameter Estimates - ' title_text], 'FontSize', 14);
    
    % Adjust plot spacing for both figures
    set(gcf, 'PaperPositionMode', 'auto');
    tight = get(gcf, 'Position');
    set(gcf, 'Position', [tight(1) tight(2) tight(3) tight(4)]);
end
