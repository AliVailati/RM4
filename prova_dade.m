%%
clc
clear all
%% es 3 e 4 
M = [87.09 9.05 0.53 0.05 0.11 0.03 0.05 0.00;
     0.48 87.32 7.72 0.46 0.05 0.06 0.02 0.02; 
     0.00 1.56 88.73 4.97 0.25 0.11 0.01 0.05; 
     0.00 0.08 3.19 86.72 3.48 0.42 0.09 0.15; 
     0.01 0.02 0.10 4.52 78.12 6.66 0.53 0.60; 
     0.00 0.02 0.06 0.15 4.54 74.73 4.81 3.18; 
     0.00 0.00 0.09 0.16 0.49 13.42 43.91 26.55;
     0.00 0.00 0.00 0.00 0.00 0.00 0.00 100.00]/100;

% Parameters
N_issuer = 100;
numRatings = size(M, 1);
cf_schedule_A = [1 0.015; 2 0.015; 3 0.015; 4 0.015; 5 1.015];
cf_schedule_BBB = [1 0.0175; 2 0.0175; 3 0.0175; 4 0.0175; 5 1.0175];
cf_schedule_BB = [1 0.025; 2 0.025; 3 0.025; 4 0.025; 5 1.025];
cf_schedule_B = [1 0.04; 2 0.04; 3 0.04; 4 0.04; 5 1.04];
cf_schedule_CCC = [1 0.06; 2 0.06; 3 0.06; 4 0.06; 5 1.06];

h = zeros(1, numRatings-1);
Prob_surv = zeros(numRatings-1, length(cf_schedule_A));
for i = 1:numRatings-1
    h(i) = -log(1 - M(i, end));
end

for i = 1:numRatings-1
    for t = 1:length(cf_schedule_A)
        Prob_surv(i, t) = exp(-h(i) * cf_schedule_A(t, 1));
    end
end
zero_rate = 0.01;
Recovery = 0.25;

%% calculate FV for the different coupons i have
FV_A = calculate_forward_values(Prob_surv, cf_schedule_A,Recovery); %Fv with coupns 1.5% sono 7 vaolri perche sono i valori di essere tra un anno in AAA, AA,....,CCC
FV_BBB = calculate_forward_values(Prob_surv, cf_schedule_BBB,Recovery);
FV_BB = calculate_forward_values(Prob_surv, cf_schedule_BB,Recovery);
FV_B = calculate_forward_values(Prob_surv, cf_schedule_B,Recovery);
FV_CCC = calculate_forward_values(Prob_surv, cf_schedule_CCC,Recovery);

disp('Forward Values:');
disp(['A: ', num2str(FV_A)]);
disp(['BBB: ', num2str(FV_BBB)]);
disp(['BB: ', num2str(FV_BB)]);
disp(['B: ', num2str(FV_B)]);
disp(['CCC: ', num2str(FV_CCC)]);

%% EV

E_FV_A = calculate_expected_forward_values(FV_A, M,3);
E_FV_BBB = calculate_expected_forward_values(FV_BBB, M,4);
E_FV_BB = calculate_expected_forward_values(FV_BB, M,5);
E_FV_B = calculate_expected_forward_values(FV_B, M,6);
E_FV_CCC = calculate_expected_forward_values(FV_CCC, M,7);

disp('Expected Forward Values:');
disp(['A: ', num2str(E_FV_A)]);
disp(['BBB: ', num2str(E_FV_BBB)]);
disp(['BB: ', num2str(E_FV_BB)]);
disp(['B: ', num2str(E_FV_B)]);
disp(['CCC: ', num2str(E_FV_CCC)]);
%% Calculate thresholds for each rating class
barriers = zeros(numRatings - 1, numRatings);
for i = 1:numRatings - 1
    cumulativeProb = 0;
    for j = 1:numRatings
        cumulativeProb = cumulativeProb + M(i, numRatings - j + 1);
        barriers(i, j) = norminv(cumulativeProb);
    end
end

disp('Barriers:');
disp(barriers);
%%
% Calculate loss matrices
Loss_A = calculate_loss(FV_A, E_FV_A, N_issuer);
Loss_BBB = calculate_loss(FV_BBB, E_FV_BBB, N_issuer);
Loss_BB = calculate_loss(FV_BB, E_FV_BB, N_issuer);
Loss_B = calculate_loss(FV_B, E_FV_B, N_issuer);
Loss_CCC = calculate_loss(FV_CCC, E_FV_CCC, N_issuer);

disp('Loss:');
disp(['A: ', num2str(Loss_A)]);
disp(['BBB: ', num2str(Loss_BBB)]);
disp(['BB: ', num2str(Loss_BB)]);
disp(['B: ', num2str(Loss_B)]);
disp(['CCC: ', num2str(Loss_CCC)]);
%% Monte Carlo Simulation of AVRs
N_sim = 1000;
% Monte Carlo Simulation of AVRs
AVR_scenario = monte_carlo_AVR(M, numRatings, N_sim, N_issuer, barriers);

%% Calculate the transition matrix from the simulation
M_simulated = zeros(numRatings - 1, numRatings);
for i = 1:numRatings - 1
    for j = 1:numRatings
        M_simulated(i, j) = mean(AVR_scenario(i, j, :)) / N_issuer;
    end
end

disp('Simulated transition matrix:');
disp(fliplr(M_simulated)); % Flip left-right for better comparison
disp('Original M matrix:');
disp(M);
%% Calculate Total Losses
Total_Loss_A = calculate_total_loss(AVR_scenario,Loss_A,numRatings,N_sim,3);
Total_Loss_BBB = calculate_total_loss(AVR_scenario,Loss_BBB,numRatings,N_sim,4);
Total_Loss_BB = calculate_total_loss(AVR_scenario,Loss_BB,numRatings,N_sim,5);
Total_Loss_B = calculate_total_loss(AVR_scenario,Loss_B,numRatings,N_sim,6);
Total_Loss_CCC = calculate_total_loss(AVR_scenario,Loss_CCC,numRatings,N_sim,7);


%% Calculate the VaR
VaR_99_A = calculate_var(Total_Loss_A);
VaR_99_BBB = calculate_var(Total_Loss_BBB);
VaR_99_BB = calculate_var(Total_Loss_BB);
VaR_99_B = calculate_var(Total_Loss_B);
VaR_99_CCC = calculate_var(Total_Loss_CCC);


disp('99% Value at Risk (VaR) for rating classes:');
disp(['A: ', num2str(VaR_99_A)]);
disp(['BBB: ', num2str(VaR_99_BBB)]);
disp(['BB: ', num2str(VaR_99_BB)]);
disp(['B: ', num2str(VaR_99_B)]);
disp(['CCC: ', num2str(VaR_99_CCC)]);

