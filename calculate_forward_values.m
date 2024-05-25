function FV = calculate_forward_values(Prob_surv, cf_schedule,Recovery)
    numRatings = size(Prob_surv, 1) + 1;
    FV = zeros(1, numRatings);
    zero_rate = 0.01;
    B = exp(-zero_rate * cf_schedule(:, 1));
    fwd_B = B(2:end) / B(1);
    
    for i = 1:numRatings-1
        FV(i) = Prob_surv(i, 2:end) * (cf_schedule(2:end, 2) .* fwd_B) + ...
                (Prob_surv(i, 1:end-1) - Prob_surv(i, 2:end)) * fwd_B * Recovery * 1;
    end
    
    FV(end) = Recovery * 1;
end
