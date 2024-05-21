function V = changeSign(V, signToConsider);
%INPUT: V matrix of the eigenvectors
% signToConsider vector with the sign to consider
%OUTPUT: V matrix of the eigenvectors with the changed sign

for i = 1:size(V, 2)
    if signToConsider(i) == 0
        V(:,i) = V(:,i);
    else 
        V(:,i) = -V(:,i);
    end
end 
end 
