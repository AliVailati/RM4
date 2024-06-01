function secondEig = second(M)

    [V_1year_sim,D_1year_sim] = eig(M);
    [D_1year_sim,ind_1year_sim] = sort(diag(D_1year_sim),'descend');
    V_1year_sim = V_1year_sim(:,ind_1year_sim);

    secondEig = V_1year_sim(:, 2); 
    if secondEig(1) <0 
        secondEig = -secondEig; %the second eigenvector is the one that simulates the trend of population
    end 
end