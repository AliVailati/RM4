function FV = calculate_forward_values2(Prob_surv, cf_schedule,Recovery, fwd_B)

    numRatings = size(Prob_surv, 1);
    FV = zeros(length(cf_schedule)-1, numRatings-2);
    
    for ii = 1:numRatings-2
        FV(:, ii) = Prob_surv(3:end, 2:end) * (cf_schedule(2:end, ii+1) .* fwd_B) + ...
                (Prob_surv(3:end, 1:end-1) - Prob_surv(3:end, 2:end)) * fwd_B * Recovery * 1;
    end
    
    FV(end, :) = Recovery * 1;
end