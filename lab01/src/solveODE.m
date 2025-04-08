function [t_sol, x_sol] = solveODE(systemODE, tspan, x0, varargin)
    % Solves a general ODE system with given parameters
    %
    % Parameters:
    %   systemODE - function handle for the ODE system (@(t,x,varargin) ...)
    %   tspan - time span for simulation [vector]
    %   x0 - initial conditions vector
    %   varargin - additional parameters to pass to the ODE function
    %
    % Returns:
    %   t_sol - time points of the solution
    %   x_sol - state vector solution
    
    % Create the ODE function which incorporates the additional parameters
    ode_func = @(t, x) systemODE(t, x, varargin{:});
    
    % Solve the ODE using ode45
    [t_sol, x_sol] = ode45(ode_func, tspan, x0);
end
