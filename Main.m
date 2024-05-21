%Final Project - RM4
%Group B - Alice Vailati, Davide Cazzaro
clear all 
close all 
clc

%% Import Data for exercise 1
addpath('Data'); 
load('M_1year.mat'); 
load('M_3year.mat'); 
load('M_5year.mat'); 

%% Remove NR
%Remove NR and distribute the value uniformly on the rows 
  M_1year = removeNR(M_1year);
  M_3year = removeNR(M_3year);
  M_5year = removeNR(M_5year);

%Add a row with all zeros except for the last element that is 1
%this row represents the absorbing state of default
M_1year = [M_1year; zeros(1,8)];
M_1year(end) = 100;
M_3year = [M_3year; zeros(1,8)];
M_3year(end) = 100;
M_5year = [M_5year; zeros(1,8)];
M_5year(end) = 100;

%% Construct the 7 years matrix: we don't have the CCC row 
M_7year = [ 38,  32.66, 6.79, 1.47, 0.35, 0.19, 0.11, 0.51, 19.93;
           1.43, 41.59, 27.29, 4.51, 0.7, 0.35, 0.03, 0.51, 23.58;
           0.06, 5.11, 48.48, 16.13, 2.28, 0.73, 0.12, 0.8, 26.28;
           0.03, 0.49, 11.01, 45.08, 7.33, 2.12, 0.34, 2.41, 31.19;
           0, 0.08, 1.2, 12.56, 24.71, 9.97, 0.94, 9.46, 41.07;
           0, 0.02, 0.26, 1.75, 8.55, 17.07, 1.84, 21.48, 49.03
         ];
M_7year = removeNR(M_7year);

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

%Last value of the seventh row taken from table 26
M_7year(7, 8) = 49.08;

diff = M_7year_simulated(7, 8)-M_7year(7, 8);

weights = M_7year(7, 1:end-1)./sum(M_7year(7, 1:end-1));
M_7year(7, 1:end-1) = M_7year(7, 1:end-1) + diff*weights;

%Normalize the matrix
M_7year = M_7year./sum(M_7year,2);
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
%%
%We want to consider the eigenvectors excluding the sign
[V_1year, V_3year, V_5year, V_7year] = signEigenvectors(V_1year, V_3year, V_5year, V_7year);

% Plot the eigenvectors 
for ii = 1:8
    plotEigenvectors(V_1year, V_3year, V_5year, V_7year, ii); 
end

%% Data exercise 2
load('M_contraction.mat'); 
load('M_expansion.mat'); 
load('issuers_contraction.mat'); 
load('issuers_expansion.mat'); 
load('M_unconditional.mat'); 
load('M_switch.mat'); 
steps_per_year = 4; % quarterly steps in a year
total_years = 5;
num_simulations = 100;
num_ratings = 8;

Me = [Me; zeros(1,8)];
Me(end) = 100;
Me = removeNR(Me);

Mc = [Mc; zeros(1,8)];
Mc(end) = 100;
Mc = removeNR(Mc);

%Normalize the matrix
Me = Me./sum(Me,2);
Mc = Mc./sum(Mc,2);

%% MC SIMULATION
M_unconditional = removeNR(M_unconditional);
%Normalize the matrix
M_unconditional = M_unconditional./sum(M_unconditional,2);

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
%% Compare the Default Probabilities
default_probs_simulated = M_avg(:, end);
default_probs_unconditional = M_unconditional(:, end);

% Plot the Default Probabilities
figure;
hold on;
plot(1:num_ratings, default_probs_simulated, 'b-o', 'LineWidth', 2, 'DisplayName', 'Simulated');
plot(1:num_ratings, default_probs_unconditional, 'r-s', 'LineWidth', 2, 'DisplayName', 'Unconditional');
title('Default Probabilities over 5 Years');
xlabel('Rating Class');
ylabel('Default Probability');
grid on;
legend('show');
hold off;

%% Exercise 3
M = [87.09 9.05 0.53 0.05 0.11 0.03 0.05 0.00;
     0.48 87.32 7.72 0.46 0.05 0.06 0.02 0.02; 
     0.00 1.56 88.73 4.97 0.25 0.11 0.01 0.05; 
     0.00 0.08 3.19 86.72 3.48 0.42 0.09 0.15; 
     0.01 0.02 0.10 4.52 78.12 6.66 0.53 0.60; 
     0.00 0.02 0.06 0.15 4.54 74.73 4.81 3.18; 
     0.00 0.00 0.09 0.16 0.49 13.42 43.91 26.55;
     0.00 0.00 0.00 0.00 0.00 0.00 0.00 100.00]/100;
numRatings = size(M, 1);
h = zeros(1, numRatings-1);
cf_schedule = [1 0.015; 2 0.015; 3 0.015; 4 0.015; 5 101.5];
zero_rate = 0.01;
B = exp(-zero_rate * cf_schedule(:, 1));
fwd_B = B(2:end) / B(1); %sono 4 
Prob_surv = zeros(numRatings-1, length(cf_schedule));
% calculate the hazard rate for every ratings
for i = 1:numRatings-1
    h(i) = -log(1 - M(i, end));
end
% calculate the survival prob for every ratings and until 5y
for i = 1:numRatings-1
    for t = 1:length(cf_schedule)
        Prob_surv(i, t) = exp(-h(i) * cf_schedule(t, 1));
    end
end
Prob_default = 1 - Prob_surv;
FV = zeros(1, numRatings);
Recovery = 0.25;

for i = 1:numRatings-1
    
   FV(i) = Prob_surv(i,2:end)*(cf_schedule(2:end,2).*fwd_B) + ...
           [Prob_surv(i,1:end-1) - Prob_surv(i,2:end)]*fwd_B*Recovery*100;
end
%Fair value conditional to be Default in 1 year
    FV(end) = Recovery * 100;
disp('Forward Values:');
disp(FV);
%now we have to make the expected values for the several FV
E_FV = zeros(numRatings - 1,1);
for i = 1:numRatings - 1
E_FV(i) = sum(FV.*M(i,:));
end
disp('Expected Forward Values:');
disp(E_FV');
