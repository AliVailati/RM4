function eigTime = eigTimeHorizon (A, B, C, D, n)
    %INPUT: 
    % A, B, C matrices of eigenvalues
    % n number of eigenvalues to be considered
    %OUTPUT: eigTime matrix with eigenvalues sorted with respect to the time horizon

    %for i from 1 to n, eig(i) is the ith eigenvalue of all the matrices
    eigTime = zeros(n,4);
    for i = 1:n
        eigTime(i,:) = [A(i), B(i), C(i), D(i)];
    end
end 