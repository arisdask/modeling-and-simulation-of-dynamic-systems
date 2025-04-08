function fig = plotInputTorque(t, inputTorqueFunc)
    % Plots the input torque u(t)
    %
    % Parameters:
    %   t - time vector
    %   inputTorqueFunc - function handle for input torque
    %
    % Returns:
    %   fig - figure handle
    
    % Calculate torque values
    u_values = zeros(size(t));
    for i = 1:1:length(t)
        u_values(i) = inputTorqueFunc(t(i));
    end
    
    fig = figure('Name', 'Input Torque');
    plot(t, u_values, 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [sec]');
    ylabel('Torque u(t) [N \cdot m]');
    title('Input Torque Function');
end
