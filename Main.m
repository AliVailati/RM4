%Final Project - RM4
%Alice Vailati, Davide Cazzaro
clear all 
close all 
clc
%% Data exercise 1
%1 year Global Corporate Average Transition Matrix
M_1year = [ 87.09, 9.05, 0.53, 0.05, 0.11, 0.03, 0.05, 0, 3.1;
            0.48, 87.32, 7.72, 0.46, 0.05, 0.06, 0.02, 0.02, 3.88;
            0.02, 1.56, 88.73, 4.97, 0.25, 0.11, 0.01, 0.05, 4.29;
            0, 0.08, 3.19, 86.72, 3.48, 0.42, 0.09, 0.15, 5.86;
            0.01, 0.02, 0.1, 4.52, 78.12, 6.66, 0.53, 0.6, 9.43;
            0, 0.02, 0.06, 0.15, 4.54, 74.73, 4.81, 3.18, 12.51;
            0, 0, 0.09, 0.16, 0.49, 13.42, 43.91, 26.55, 15.39
           ]; 

M_3year = [ 65.54, 22.15, 2.32, 0.32, 0.26, 0.08, 0.11, 0.13, 9.08;
            1.11, 67.26, 18.04, 1.92, 0.32, 0.2, 0.03, 0.11, 11.01;
            0.05, 3.67, 70.68, 11.14, 1.1, 0.38, 0.08, 0.22, 12.67;
            0.02, 0.24, 7.9, 66.78, 6.71, 1.42, 0.25, 0.74, 15.93;
            0.01, 0.05, 0.43, 10.32, 49.13, 11.3, 1.15, 3.39, 24.23;
            0, 0.02, 0.16, 0.63, 9.08, 42.39, 5.26, 11.56, 30.9;
            0, 0, 0.11, 0.5, 1.51, 16.52, 9.73, 42.29, 29.35
           ];


M_5year = [ 49.52, 28.83, 4.75, 0.8, 0.34, 0.16, 0.08, 0.34, 15.17;
            1.37, 52.56, 23.97, 3.38, 0.53, 0.35, 0.04,  0.3, 17.49;
            0.06, 4.7, 57.8, 14.39, 1.8, 0.6, 0.13, 0.46, 20.05;
            0.02, 0.37, 10.14, 53.82, 7.5, 1.94, 0.34, 1.58,  24.29;
            0.01, 0.07, 0.86, 12.27, 33.54, 11.15, 1.15, 6.51, 34.45;
            0.01, 0.02, 0.2, 1.26, 9.59, 25.75,  3.21, 17.4, 42.56;
            0, 0, 0.09, 0.68, 2.49, 12.34, 2.7, 46.35,  35.35
           ];

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
M_1year(end) = 1;
M_3year = [M_3year; zeros(1,8)];
M_3year(end) = 1;
M_5year = [M_5year; zeros(1,8)];
M_5year(end) = 1;

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

% Plot the second eigenvectors
figure
plot(V_1year(1:end-1,2),'-o','DisplayName','1 year')
 hold on
plot(V_3year(:,2),'-square','DisplayName','3 year')
plot(V_5year(:,2),'-diamond','DisplayName','5 year')
%% Data exercise 2
Me =  [ 98.21, 1.66, 0.11, 0.02, 0.02, 0, 0, 0;
        0.15, 98.08, 1.61, 0.12, 0.01, 0.03, 0.01, 0;
        0.02, 0.53, 98.06, 1.21, 0.11, 0.06, 0, 0;
        0.01, 0.07, 1.47, 96.94, 1.25, 0.22, 0.02, 0.02;
        0.01, 0.03, 0.19, 1.93, 95.31, 2.25, 0.16, 0.12;
        0, 0.02, 0.07, 0.1, 1.7, 95.91, 1.31, 0.88;
        0.05, 0, 0.19, 0.23, 0.47, 3.57, 87.32, 8.17
      ];
Me_issuers = [6581, 19458, 36404, 24529, 18161, 20002, 2129];

Mc =  [ 97.99, 1.76, 0.25, 0, 0, 0, 0, 0;
        0.18, 96.89, 2.79, 0.05, 0.09, 0, 0, 0;
        0.02, 0.88, 96.44, 2.59, 0.07, 0, 0, 0;
        0.04, 0.04, 1.11, 96.31, 2.33, 0.07, 0, 0.11;
        0, 0.06, 0.06, 1.39, 94.98, 2.72, 0.42, 0.36;
        0, 0.06, 0.06, 0.11, 0.72, 95.02, 2.27, 1.77;
        0, 0, 0, 0, 0, 1.2, 85.6, 13.2
      ];
Mc_issuers = [795, 2186, 4330, 2708, 1655, 1806, 250];