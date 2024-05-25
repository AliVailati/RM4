function Loss = calculate_loss(FV, E_FV, N_issuer)
    numRatings = length(FV);
    Loss = zeros(1, numRatings);

        for j = 1:numRatings
            Loss(j) = (FV(numRatings - j + 1) - E_FV) / N_issuer;
        end
end
