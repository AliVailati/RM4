function eigTime = eigTimeHorizon (A, B, C, D, n)
    %INPUT: 
    %        A, B, C, D = matrices of eigenvalues
    %        n = number of eigenvalues to be considered
    %OUTPUT: 
    %        eigTime = matrix with eigenvalues sorted with respect to the time horizon

    %The objective of this function is to cluster the eigenvalues of the
    %different matrices so that we have the i-th eigenvalue of each marix
    %in the same row of eigTime

    eigTime = zeros(n,4);
    for i = 1:n
        eigTime(i,:) = [A(i), B(i), C(i), D(i)];
    end
end 