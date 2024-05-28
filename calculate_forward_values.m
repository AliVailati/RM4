function [FV, E_FV] = calculate_forward_values(Prob_surv, cf_schedule,Recovery, fwd_B, Notional, M)
    %Prob_surv = matrix of survival probabilities: 7x5 matrix
    %cf_schedule = 5x6 but with the first column with the years
    %recovery = scalar
    %fwd_B = 4x1
    numRatings = size(Prob_surv, 1);
    numRateCoupons = size(cf_schedule, 1); 
    FV = zeros(numRatings, numRateCoupons);
    E_FV = zeros(1, numRateCoupons);
    for k = 1:numRateCoupons
        FV(:, k) = Prob_surv(:, 2:end) * (cf_schedule(2:end, k+1) .* fwd_B) + ...
                (Prob_surv(:, 1:end-1) - Prob_surv(:, 2:end)) * fwd_B * Recovery * Notional;
    end
    FV = FV'; 
    FV = [FV, zeros(numRateCoupons, 1)]; 
    FV(:, end) = Recovery*Notional;

   for k = 1:numRateCoupons
       E_FV(k) = sum(FV(k, :).*M(k+2, :));
   end
end