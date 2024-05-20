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

%Nella 7 anni manca la riga della CCC
%M_7year = [ 38,  32.66, 6.79, 1.47, 0.35, 0.19, 0.11, 0.51, 19.93;
%            1.43, 41.59, 27.29, 4.51, 0.7, 0.35, 0.03, 0.51, 23.58;
%            0.06, 5.11, 48.48, 16.13, 2.28, 0.73, 0.12, 0.8, 26.28;
%            0.03, 0.49, 11.01, 45.08, 7.33, 2.12, 0.34, 2.41, 31.19;
%            0, 0.08, 1.2, 12.56, 24.71, 9.97, 0.94, 9.46, 41.07;
%            0, 0.02, 0.26, 1.75, 8.55, 17.07, 1.84, 21.48, 49.03
%          ];

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

%Normalize the matrix
M_1year = M_1year./sum(M_1year,2);
M_3year = M_3year./sum(M_3year,2);
M_5year = M_5year./sum(M_5year,2);
%%
%Analysis of the Eigenvalues of the transition matrix
%1 year
[V_1year,D_1year] = eig(M_1year);
%3 year
[V_3year,D_3year] = eig(M_3year);
%5 year
[V_5year,D_5year] = eig(M_5year);

%We have to sort the eigenvalues and the eigenvectors
[D_1year,ind_1year] = sort(diag(D_1year),'descend');
V_1year = V_1year(:,ind_1year);

[D_3year,ind_3year] = sort(diag(D_3year),'descend');
V_3year = V_3year(:,ind_3year);

[D_5year,ind_5year] = sort(diag(D_5year),'descend');
V_5year = V_5year(:,ind_5year);

%extract, grouped by years, the first n eigenvalues for each matrix
eigTime = eigTimeHorizon (D_1year, D_3year, D_5year, 4);

%Plot the natural logarithm of eigenvalues with respect to the time horizon
figure
plot([1,3,5],log(eigTime(2, :)),'-o','DisplayName','Second Eigenvalue')
hold on
plot([1,3,5],log(eigTime(3, :)),'-square','DisplayName','Third Eigenvalue')
plot([1,3,5],log(eigTime(4, :)),'-diamond','DisplayName','Fourth Eigenvalue')
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
