function state_paths = MCpath (steps_per_year, total_years, num_simulations, M_switch)

total_steps = steps_per_year * total_years;
state_paths = zeros(num_simulations, total_steps);

rng(0);

% Monte Carlo Simulation of State Paths
for sim = 1:num_simulations
    %randomly starting in expansion or recession
    state_paths(sim, :) = createPath(total_steps, M_switch);
end