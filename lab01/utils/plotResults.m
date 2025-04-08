function plotResults(t_true, x_true, t_est, x_est, e_q, e_q_dot, title_str)
    % Creates a 2x2 figure with trajectory comparisons and errors
    %
    % Inputs:
    %   t_true - True time vector
    %   x_true - True state trajectories
    %   t_est - Time vector for the estimated trajectory
    %   x_est - Estimated state trajectories
    %   e_q - Error in angle
    %   e_q_dot - Error in angular velocity
    %   title_str - Title for the figure
    
    figure;
    subplot(2,2,1);
    plot(t_true, x_true(:,1), 'b-', t_est, x_est(:,1), 'r--', 'LineWidth', 1.5);
    title('Angle q(t)');
    xlabel('Time [sec]');
    ylabel('Angle [rad]');
    legend('True', 'Estimated');
    grid on;
    
    subplot(2,2,2);
    plot(t_true, x_true(:,2), 'b-', t_est, x_est(:,2), 'r--', 'LineWidth', 1.5);
    title("Angular Velocity q'(t)");
    xlabel('Time [sec]');
    ylabel('Angular Velocity [rad/s]');
    legend('True', 'Estimated');
    grid on;
    
    subplot(2,2,3);
    plot(t_est, e_q, 'k-', 'LineWidth', 1.5);
    title('Error in Angle e_q(t)');
    xlabel('Time [sec]');
    ylabel('Error [rad]');
    grid on;
    
    subplot(2,2,4);
    plot(t_est, e_q_dot, 'k-', 'LineWidth', 1.5);
    title("Error in Angular Velocity e_q'(t)");
    xlabel('Time [sec]');
    ylabel('Error [rad/s]');
    grid on;
    
    sgtitle(title_str);
end
