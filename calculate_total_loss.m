function total_loss = calculate_total_loss(AVR_scenario, Loss, numRatings, N_sim, flag)
    total_loss = zeros(N_sim, 1);
    for j = 1:N_sim
        total_loss_scenario = 0;
        for k = 1:numRatings
            total_loss_scenario = total_loss_scenario + AVR_scenario(flag, k, j) * Loss(k);
        end
        total_loss(j) = total_loss_scenario;
    end
end
