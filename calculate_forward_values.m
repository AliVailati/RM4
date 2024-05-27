function FV = calculate_forward_values(Prob_surv, cf_schedule,Recovery, fwd_B)

    numCoupons = length(cf_schedule)-1;
    numRatings = numCoupons;
    FV = zeros(numCoupons, numRatings);
    
    for ii = 1 : numCoupons
        FV(ii, :) = Prob_surv(3:end, 2:end) * (cf_schedule(2:end, ii+1) .* fwd_B) + ...
                (Prob_surv(3:end, 1:end-1) - Prob_surv(3:end, 2:end)) * fwd_B * Recovery * 1;
    end
    FV = [FV, zeros(numRatings, 1)]; 
    FV(:, end) = Recovery * 1;
end