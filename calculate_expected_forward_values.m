function E_FV = calculate_expected_forward_values(FV, M, flag)
    E_FV = sum(FV .* M(flag, :));
end
