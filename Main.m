%Final Project - RM4
%Group B - Alice Vailati, Davide Cazzaro
clear all 
close all 
clc
tic

%% EXERCISE 1 -> Import Data and functions
addpath('Data'); 
load('M_1year.mat'); 
load('M_3year.mat'); 
load('M_5year.mat'); 
load('M_7year.mat'); 

addpath('Function_ex1'); 

%% Remove NR
%Remove NR and distribute the value uniformly on the rows 
  M_1year = removeNR(M_1year);
  M_3year = removeNR(M_3year);
  M_5year = removeNR(M_5year);
  M_7year = removeNR(M_7year);

%Add a row with all zeros except for the last element that is 1
%this row represents the absorbing state of default
M_1year = [M_1year; zeros(1,8)];
M_1year(end) = 100;
M_3year = [M_3year; zeros(1,8)];
M_3year(end) = 100;
M_5year = [M_5year; zeros(1,8)];
M_5year(end) = 100;

%% Construction of the 7 years matrix -> there is not the CCC row 
%Add a row of zeros in which we will add values later for row CCC
M_7year = [M_7year; zeros(1,8)];

%Add a row of zeros except for the default case equal to 1 (default is an absorbing state)
M_7year = [M_7year; zeros(1,8)];
M_7year(end) = 100;

%Normalize the matrix
M_1year = M_1year./sum(M_1year,2);
M_3year = M_3year./sum(M_3year,2);
M_5year = M_5year./sum(M_5year,2);

%Simulate the 7 years matrix and take it back to have the sum of the rows equal to 100
M_7year_simulated = 100*M_5year*M_1year^2;

%Copy the values of the seventh row of the simulated matrix in the 7 years matrix
M_7year(end-1,:) = M_7year_simulated(end-1,:);

%Last value of the seventh row is taken from table 26
M_7year(7, 8) = 49.08;

diff = M_7year_simulated(7, 8)-M_7year(7, 8);

weights = M_7year(7, 1:end-1)./sum(M_7year(7, 1:end-1));
M_7year(7, 1:end-1) = M_7year(7, 1:end-1) + diff*weights;

%Normalize the matrix
M_7year = M_7year./sum(M_7year,2);

clear diff weights
%%
%Analysis of the Eigenvalues of the transition matrix
%1 year
[V_1year,D_1year] = eig(M_1year);
%3 year
[V_3year,D_3year] = eig(M_3year);
%5 year
[V_5year,D_5year] = eig(M_5year);
%7 year
[V_7year,D_7year] = eig(M_7year);

%% Eigenvalues
%We have to sort the eigenvalues and the eigenvectors
[D_1year,ind_1year] = sort(diag(D_1year),'descend');
V_1year = V_1year(:,ind_1year);

[D_3year,ind_3year] = sort(diag(D_3year),'descend');
V_3year = V_3year(:,ind_3year);

[D_5year,ind_5year] = sort(diag(D_5year),'descend');
V_5year = V_5year(:,ind_5year);

[D_7year,ind_7year] = sort(diag(D_7year),'descend');
V_7year = V_7year(:,ind_7year);

%extract, grouped by years, the first n eigenvalues for each matrix
eigTime = eigTimeHorizon (D_1year, D_3year, D_5year, D_7year, 4);

%Plot the natural logarithm of eigenvalues with respect to the time horizon
figure
plot([1,3,5,7],log(eigTime(2, :)),'-o','DisplayName','Second Eigenvalue')
hold on
plot([1,3,5,7],log(eigTime(3, :)),'-square','DisplayName','Third Eigenvalue')
plot([1,3,5,7],log(eigTime(4, :)),'-diamond','DisplayName','Fourth Eigenvalue')
xlabel('Time Horizon')
ylabel('Logarithm of Eigenvalues')
title('Eigenvalues with respect to the Time Horizon')
legend('show', 'Location', 'best')
grid on

clear ind_1year ind_3year ind_5year ind_7year
%% Change sign to eigenvectors
[V_1year, V_3year, V_5year, V_7year] = signEigenvectors(V_1year, V_3year, V_5year, V_7year);

% Plot the eigenvectors 
for ii = 1:8
    plotEigenvectors(V_1year, V_3year, V_5year, V_7year, ii); 
end

clear ii
%% EXERCISE 2 -> Import Data and functions
addpath('Function_ex2');
load('Mc.mat'); 
load('Me.mat'); 
load('Mc_issuers.mat'); 
load('Me_issuers.mat'); 
load('M_unconditional.mat'); 
load('M_switch.mat'); 
steps_per_year = 4; % quarterly steps in a year
total_years = 5;
num_simulations = 10000;
num_ratings = 8;

%% Add row for default and remove NR 
%In this case removeNR is useful to have all the row summing up to 100
Me = [Me; zeros(1,8)];
Me(end) = 100;
Me = removeNR(Me);

Mc = [Mc; zeros(1,8)];
Mc(end) = 100;
Mc = removeNR(Mc);

M_unconditional = removeNR(M_unconditional);

%Normalize the matrices
Me = Me./sum(Me,2);
Mc = Mc./sum(Mc,2);
M_unconditional = M_unconditional./sum(M_unconditional,2);

%% MC SIMULATION
total_steps = steps_per_year * total_years;
state_paths = zeros(num_simulations, total_steps);

rng(0);

% Monte Carlo Simulation of State Paths
for sim = 1:num_simulations
    %randomly starting in expansion or recession
    state_paths(sim, :) = createPath(total_steps, M_switch);
end

%% Calculate the 5-year Transition Matrices -> 3d matrix
transition_matrices = zeros(num_ratings, num_ratings, num_simulations); % I need a 3d matrix -> 3rd dim the simulations

% Compute the Transition Matrices for each Simulation
for sim = 1:num_simulations
    transition_matrices(:, :, sim) = computeTransitionMatrix(state_paths(sim, :), num_ratings, Me, Mc);
end

clear sim
%% Compute the Average 5-year Transition Matrix
M_avg = mean(transition_matrices, 3);

disp('Average 5-year Simulated Transition Matrix:');
disp(M_avg);

%% MSE between simulated and unconditional 
M_diff = M_unconditional-M_avg; 
M_diff = M_diff.^2; 
MSE = mean(M_diff(:)); 
disp ('The MSE between the unconditional one and the simulated is:'); 
disp(MSE); 

clear M_diff
%% Compare the Default Probabilities
default_probs_simulated = M_avg(:, end);
default_probs_unconditional = M_unconditional(:, end);

% Plot the Default Probabilities
figure;
hold on;
plot(1:num_ratings, default_probs_simulated, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
plot(1:num_ratings, default_probs_unconditional, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Unconditional');
title('Default Probabilities over 5 Years');
xlabel('Rating Class');
ylabel('Default Probability');
grid on;
legend('show', 'Location', 'best');
hold off;

%% Exercise 3 -> Import data and functions
addpath('Function_ex3');
load('coupon_rates.mat');
%Consider the one year unconditional matrix
M = M_1year;
%Annual payment from year 1 to 5 
coupon_times = [1 2 3 4 5]; 
%coupons for the initial ratings 
zero_rate = 0.01;
Recovery = 0.25; 
N_issuers = 100; 
Notional = 1; 
%% Consider all the 5 portfolios -> 5 issuers, cf_schedule is a matrix 5x6
cf_schedule = coupons(coupon_times, coupon_rates);

%calculate the discount factor for each year
B = exp(-zero_rate * cf_schedule(:, 1));
fwd_B = B(2:end) / B(1);

%hazard rate 
h = -log(1 - M(1:end-1, end));

% calculate the survival prob for every ratings and until 5y
probSurv = exp(-h*coupon_times);
probDefault = 1 - probSurv;
%%
%calculate the forward value 
[FV, E_FV] = calculate_forward_values(probSurv, cf_schedule ,Recovery, fwd_B, Notional, M);

% Define the column width for alignment
columnWidth = 8; 

% Display the forward values
fprintf('Forward Values:\n');
% Print the header row with years
header = ['     ', repmat(['%', num2str(columnWidth), 's '], 1, size(FV, 2))];
fprintf(header, 'AAA', 'AA', 'A', 'BBB', 'BB', 'B', 'CCC', 'D');
fprintf('\n');
fprintf(['A:   ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(FV, 2)), '\n'], FV(1, :));
fprintf(['BBB: ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(FV, 2)), '\n'], FV(2, :));
fprintf(['BB:  ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(FV, 2)), '\n'], FV(3, :));
fprintf(['B:   ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(FV, 2)), '\n'], FV(4, :));
fprintf(['CCC: ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(FV, 2)), '\n'], FV(5, :));
fprintf('\n');

% Define the column width for alignment
columnWidth = 8; 

% Display the forward values
fprintf('Expected Forward Values:\n');
% Print the header row with years
fprintf('\n');
fprintf(['A:   ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(E_FV, 1)), '\n'], E_FV(1));
fprintf(['BBB: ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(E_FV, 1)), '\n'], E_FV(2));
fprintf(['BB:  ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(E_FV, 1)), '\n'], E_FV(3));
fprintf(['B:   ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(E_FV, 1)), '\n'], E_FV(4));
fprintf(['CCC: ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(E_FV, 1)), '\n'], E_FV(5));

%questi sono i PV at today di un bond AAA,AA,A,ecc...
%% Losses
E_FV = E_FV'; 
E_FV = E_FV*ones(1, 8); 
Loss = -(FV-E_FV)/N_issuers;

columnWidth = 8; 

% Display the forward values
fprintf('Losses:\n');
% Print the header row with years
fprintf('\n');
fprintf(['A:   ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(Loss, 2)), '\n'], Loss(1, :));
fprintf(['BBB: ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(Loss, 2)), '\n'], Loss(2, :));
fprintf(['BB:  ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(Loss, 2)), '\n'], Loss(3, :));
fprintf(['B:   ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(Loss, 2)), '\n'], Loss(4, :));
fprintf(['CCC: ', repmat(['%', num2str(columnWidth), '.4f '], 1, size(Loss, 2)), '\n'], Loss(5, :));
% % calculate the treshold for every rating class
% barriers = zeros(numRatings - 1, numRatings);
% cumulativeProb = 0;
% for i = 1:numRatings - 1 %righe
%     cumulativeProb = 0;
%     for j = 1:numRatings % colonne
%         cumulativeProb = cumulativeProb + M(i,numRatings - j + 1); %devo andare al contrario -> vedere grafico avr
%         barriers(i,j) = norminv(cumulativeProb);
%     end
% end 
% disp('Barriers:');
% disp(barriers);
% 
% %calculation of the Loss for every case
% Loss_matrix = zeros(numRatings - 1,numRatings);
% for i = 1:numRatings - 1
%     for j = 1:numRatings
%         Loss_matrix(i,j) = (FV(numRatings - j + 1) - E_FV(i)) / N_issuer; %in verita dovrebbero essere i returns
%     end
% end
% disp('Loss:');
% disp(Loss_matrix);
% %prima  se vado in def via via a crescere
% 
% 
% 
% 
% 
% %% Monte Carlo Simulation of AVRs
% N_sim = 1000;
% rho = zeros(1, 8);
% for i = 1:8
%     rho(i) = rho_R(M, i);
% end
% %devo fare una matrice 3d per poter tenere conto di tutte le varie
% %simulazioni
% AVR_scenario = zeros(numRatings - 1, numRatings, N_sim);
% 
% for i = 1:numRatings - 1 %riempio le righe
%     Y = randn(N_sim, 1); % Simulo il common factor una sola volta per rating class; giusto?
%     for j = 1:N_sim %tridimensionalità
%         for i = 1:N_issuer
%             epsilon = randn; % simulo ogni volta -> è rischio idiosincratico
%             V = rho(i) * Y(j) + sqrt(1 - rho(i)^2) * epsilon;
%             for k = 1:numRatings %riempio le colonne
%                 if k == 1
%                     if V <= barriers(i, k) %guardo default
%                         AVR_scenario(i, k, j) = AVR_scenario(i, k, j) + 1;
%                         break;
%                     end
%                 else
%                     if V <= barriers(i, k) && V > barriers(i, k-1)
%                         AVR_scenario(i, k, j) = AVR_scenario(i, k, j) + 1;
%                         break;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% %% Calculate the transition matrix from the simulation
% M_simulated = zeros(numRatings - 1, numRatings);
% for i = 1:numRatings - 1
%     for j = 1:numRatings
%         M_simulated(i, j) = mean(AVR_scenario(i, j, :)) / N_issuer;
%     end
% end
% 
% disp('Simulated transition matrix:');
% disp(fliplr(M_simulated)); %inverto l'ordine degli elementi di ogni riga cosi confronto meglio
% disp('Original M matrix:');
% disp(M);
% 
% %% Calculate VaR at 99%
% total_losses = zeros(N_sim, numRatings - 1);%AAA AA A BBB BB B CCC 
% for i = 1:numRatings - 1 %colonne
%     for j = 1:N_sim
%         total_loss = 0;
%         for k = 1:numRatings %serve per usare AVR_scenario
%             total_loss = total_loss + AVR_scenario(i, k, j) * Loss_matrix(i, k);
%         end
%         total_losses(j, i) = total_loss;
%     end
% end
% 
% VaR_99 = zeros(numRatings - 1, 1);
% for i = 1:numRatings - 1
%     VaR_99(i) = quantile(total_losses(:, i), 0.99);
% end
% 
% disp('99% Value at Risk (VaR) for each rating class:');
% disp(VaR_99);
% 
% %% Plot the PDF of total losses for each rating class
% % Define the rating class names
% ratingClassNames = {'AAA', 'AA', 'A', 'BBB', 'BB', 'B', 'CCC'};
% 
% % Plot the PDF of total losses for each rating class
% figure;
% hold on;
% colors = lines(numRatings - 1);
% for r = 1:numRatings - 1
%     [f, xi] = ksdensity(total_losses(:, r));
%     plot(xi, f, 'Color', colors(r, :), 'DisplayName', ratingClassNames{r});
% end
% title('PDF of Total Losses for Different Rating Classes');
% xlabel('Total Loss');
% ylabel('Density');
% legend show;
% hold off;
% 
% %interessante vedere come cambia rispetto alle varie rating class!
% toc 