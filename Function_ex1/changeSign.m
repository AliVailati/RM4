function V = changeSign(V, signToConsider)
    %INPUT: 
    %       V = matrix of the eigenvectors
    %       signToConsider = vector with the sign to consider
    %OUTPUT: 
    %       V = matrix of the eigenvectors with the changed sign

    %The objective of this function is to change the sign of the
    %eigenvectors if needed

for i = 1:size(V, 2)
    if signToConsider(i) == 0
        V(:,i) = V(:,i);
    else 
        V(:,i) = -V(:,i);
    end
end 
end 
