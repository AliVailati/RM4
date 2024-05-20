%Final Project - RM4
%Alice Vailati, Davide Cazzaro
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

%Add a row of zeros in which we will add values later
M_7year = [M_7year; zeros(1,8)];

%Add a row of zeros except for the default case equal to 1
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

% Plot the eigenvectors
figure 
subplot(4, 2, 1)
plot(V_1year(:,1),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,1),'-square','DisplayName','3 year')
plot(V_5year(:,1),'-diamond','DisplayName','5 year')
grid on
title ('First eigenvector')
hold off

subplot(4, 2, 2)
plot(V_1year(:,2),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,2),'-square','DisplayName','3 year')
plot(V_5year(:,2),'-diamond','DisplayName','5 year')
grid on
title ('Second eigenvector')
hold off

subplot(4, 2, 3)
plot(V_1year(:,3),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,3),'-square','DisplayName','3 year')
plot(V_5year(:,3),'-diamond','DisplayName','5 year')
grid on
title ('Third eigenvector')
hold off

subplot(4, 2, 4)
plot(V_1year(:,4),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,4),'-square','DisplayName','3 year')
plot(V_5year(:,4),'-diamond','DisplayName','5 year')
grid on
title ('Fourth eigenvector')
hold off

subplot(4, 2, 5)
plot(V_1year(:,5),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,5),'-square','DisplayName','3 year')
plot(V_5year(:,5),'-diamond','DisplayName','5 year')
grid on
title ('Fifth eigenvector')
hold off

subplot(4, 2, 6)
plot(V_1year(:,6),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,6),'-square','DisplayName','3 year')
plot(V_5year(:,6),'-diamond','DisplayName','5 year')
grid on
title ('Sixth eigenvector')
hold off

subplot(4, 2, 7)
plot(V_1year(:,7),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,7),'-square','DisplayName','3 year')
plot(V_5year(:,7),'-diamond','DisplayName','5 year')
grid on
title ('Seventh eigenvector')
hold off

subplot(4, 2, 8)
plot(V_1year(:,8),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,8),'-square','DisplayName','3 year')
plot(V_5year(:,8),'-diamond','DisplayName','5 year')
grid on
title ('Eigth eigenvector')
hold off

%% Data exercise 2
load('M_contraction.mat'); 
load('M_expansion.mat'); 
load('issuers_contraction.mat'); 
load('issuers_expansion.mat'); 

%% es 2 MC SIMULATION
M_unconditional = [49.52 28.83 4.75 0.80 0.34 0.16 0.08 0.34 15.17
                    1.37 52.56 23.97 3.38 0.53 0.35 0.04 0.30 17.49
                    0.06 4.70 57.80 14.39 1.80 0.60 0.13 0.46 20.05
                    0.02 0.37 10.14 53.82 7.50 1.94 0.34 1.58 24.29
                    0.01 0.07 0.86 12.27 33.54 11.15 1.15 6.51 34.45
                    0.01 0.02 0.20 1.26 9.59 25.75 3.21 17.40 42.56
                    0.00 0.00 0.09 0.68 2.49 12.34 2.70 46.35 35.35
                    0 0 0 0 0 0 0 100 0 ]/100;
M_unconditional = removeNR(M_unconditional);
%regime switching matrix
M_switch = [0.85, 0.15;
            0.692, 0.308];

steps_per_year = 4; % quarterly steps in a year
total_years = 5;
total_steps = steps_per_year * total_years;
num_simulations = 100;
num_ratings = 8;
state_paths = zeros(num_simulations, total_steps);

rng(0);

% Monte Carlo Simulation of State Paths
for sim = 1:num_simulations
    %randomly starting in expansion or recession
    current_state = randi([0, 1]); %expansion (0) or recession (1)
    for step = 1:total_steps
        if current_state == 0
            if rand < M_switch(1, 1) %rand gives me int betw [0,1]
                next_state = 0; % stay in expansion
            else
                next_state = 1; % switch to recession
            end
        else
            if rand < M_switch(2, 1)
                next_state = 0; % switch to expansion
            else
                next_state = 1; % stay in recession
            end
        end
        state_paths(sim, step) = next_state;
        current_state = next_state;
    end
end

%% Calculate the 5-year Transition Matrices -> 3d matrix
transition_matrices = zeros(num_ratings, num_ratings, num_simulations); % i need a 3d matrix -> 3rd dim the simulations
for sim = 1:num_simulations
    transition_matrix = eye(num_ratings);
    for step = 1:total_steps
        if state_paths(sim, step) == 0
            transition_matrix = transition_matrix * Me;
        else
            transition_matrix = transition_matrix * Mc;
        end
    end
    transition_matrices(:, :, sim) = transition_matrix;
end

%% Compute the Average 5-year Transition Matrix
avg_transition_matrix = mean(transition_matrices, 3);
disp('Average 5-year Simulated Transition Matrix:');
disp(avg_transition_matrix);

%% Plot Results
default_probs_simulated = avg_transition_matrix(:, end);
default_probs_unconditional = M_unconditional(:, end);

% def prob
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
