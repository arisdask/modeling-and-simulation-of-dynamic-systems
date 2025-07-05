function dfilter_states_dt = updateFilterStates(filter_states, x_true, u, lambda)
    % Update all filter states
    % Each filter implements: G(s) = 1/(s + λ)
    % State space: dx_f/dt = -λ*x_f + input, output = x_f
    
    dfilter_states_dt = zeros(3, 1);
    
    % Filter for x1
    dfilter_states_dt(1) = -lambda * filter_states(1) + x_true(1);
    
    % Filter for x2
    dfilter_states_dt(2) = -lambda * filter_states(2) + x_true(2);
    
    % Filter for input u
    dfilter_states_dt(3) = -lambda * filter_states(3) + u;
end
