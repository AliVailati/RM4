function [V_1year, V_3year, V_5year, V_7year, signToConsider] = signEigenvectors(V_1year, V_3year, V_5year, V_7year)
    %INPUT: 
    %        V_1year, V_3year, V_5year, V_7year = matrix of the eigenvectors
    %OUTPUT: 
    %        V_1year, V_3year, V_5year, V_7year = matrix of the
    %        eigenvectors with changed sign

    %The objective of this function is to made the eigenvectors comparable
    %making all of them to start from a positive value

%Sign of the first row of the eigenvectors
signEigenvectors = [sign(V_1year(1,:)); sign(V_3year(1,:)); sign(V_5year(1,:)); sign(V_7year(1,:))];

pos = ones(1, 8); 
signPositive = sign(pos);
%Subtract to all of this the value of 1 in order to obtain a matrix of
%zeros (same sign as 1) and -1 (different sign from 1)

signToConsider = signEigenvectors-ones(4,1)*signPositive;
signToConsider = sign(signToConsider);

%Change the sign of the eigenvectors
V_1year = changeSign(V_1year, signToConsider(1,:));
V_3year = changeSign(V_3year, signToConsider(2,:));
V_5year = changeSign(V_5year, signToConsider(3,:));
V_7year = changeSign(V_7year, signToConsider(4,:));

end