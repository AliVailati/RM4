function barriers = calculate_barriers(M)
    numRatings = size(M, 1);
    barriers = zeros(numRatings - 1, numRatings);
    
    for i = 1:numRatings - 1
        cumulativeProb = 0;
        for j = 1:numRatings
            cumulativeProb = cumulativeProb + M(i, numRatings - j + 1);
            barriers(i, j) = norminv(cumulativeProb);
        end
    end
end
