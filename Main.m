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
M_7year = [ 38,  32.66, 6.79, 1.47, 0.35, 0.19, 0.11, 0.51, 19.93;
            1.43, 41.59, 27.29, 4.51, 0.7, 0.35, 0.03, 0.51, 23.58;
            0.06, 5.11, 48.48, 16.13, 2.28, 0.73, 0.12, 0.8, 26.28;
            0.03, 0.49, 11.01, 45.08, 7.33, 2.12, 0.34, 2.41, 31.19;
            0, 0.08, 1.2, 12.56, 24.71, 9.97, 0.94, 9.46, 41.07;
            0, 0.02, 0.26, 1.75, 8.55, 17.07, 1.84, 21.48, 49.03
          ];
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