function state_paths = createPath(total_steps, M_switch)
    %INPUT : 
    %total_steps = number of steps that we have to consider
    %M_switch = matrix of switching state 
    %OUTPUT = vector of path

    state_paths = zeros(total_steps, 1); 
    state_paths(1) = randi([0, 1]); %expansion (0) or recession (1)
    
    for step = 1:total_steps-1
        currentState = state_paths(step); 
        if currentState == 0
            switchProb = M_switch(1, 2);  %prob of switching from expansion to recession
        else 
            switchProb = M_switch(2, 1);  %prob of switching from recession to expansion
        end 
        if rand < switchProb   
            state_paths(step+1) = 1 - currentState; % switch from current state to the other state
        else
            state_paths(step+1) = currentState; % stay in the current state
        end
    end