function T =  computeTransitionMatrix(state_path, num_ratings, Me, Mc)
    %INPUT
    %       state_path = a vector of states
    %       num_ratings = number of ratings
    %       Me = matrix of expansion
    %       Mc = matrix of contraction
    %OUTPUT
    % T: transition matrix

    %The objective of this function is to compute the Transition Matrix for
    %each path multiplying an eye matrix for Me each time in which thebpath
    %is in expansion and for Mc when is in contraction
    
    T = eye(num_ratings, num_ratings);
    for i = 1:length(state_path)-1
        if state_path(i) == 0
            T = T * Me;
        else
            T = T * Mc;
        end
    end
end 