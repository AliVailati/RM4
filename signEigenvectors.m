function [V_1year, V_3year, V_5year, V_7year] = signEigenvectors(V_1year, V_3year, V_5year, V_7year);
%INPUT: V_1year, V_3year, V_5year, V_7year matrix of the eigenvectors
%OUTPUT: V_1year, V_3year, V_5year, V_7year matrix of the eigenvectors with the same sign

%Sign of the first row of the eigenvectors
signEigenvectors = [sign(V_1year(1,:)); sign(V_3year(1,:)); sign(V_5year(1,:)); sign(V_7year(1,:))];

%Consider the fisrt one as a reference and subtract it to all the others
signToConsider = signEigenvectors-ones(4,1)*signEigenvectors(1,:);
signToConsider = sign(signToConsider);

%Change the sign of the eigenvectors
V_3year = changeSign(V_3year, signToConsider(2,:));
V_5year = changeSign(V_5year, signToConsider(3,:));
V_7year = changeSign(V_7year, signToConsider(4,:));

end