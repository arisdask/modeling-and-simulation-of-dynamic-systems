function fig = plotPendulumStates(t, x)
    % Plots the pendulum state 
    % [angle q(t) = x1, angular velocity q'(t) = x2]
    %
    % Parameters:
    %   t - time vector
    %   x - state vector [q(t), q'(t)]
    %
    % Returns:
    %   fig - figure handle
    
    fig = figure('Name', 'Pendulum States');
    
    % Plot angle q(t)
    subplot(2, 1, 1);
    plot(t, x(:,1), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [sec]');
    ylabel('Angle q(t) [rad]');
    title('Pendulum Angle q(t)');
    
    % Plot angular velocity q'(t)
    subplot(2, 1, 2);
    plot(t, x(:,2), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [sec]');
    ylabel('Angular Velocity [rad/sec]');
    title("Pendulum Angular Velocity q'(t)");
end
