function VaR_99 = calculate_var(total_loss)
    VaR_99 = quantile(total_loss, 0.99);
end
