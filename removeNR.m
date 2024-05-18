function M = removeNR(M_old)
    %INPUT
    % M_old: matrix with NR
    %OUTPUT
    % M: matrix without NR

    %NR is the last column of M_old
    NR = M_old(:,end);
    %Remove NR from M_old
    M_old = M_old(:,1:end-1);
    %Distribute NR uniformly on the row
    M = M_old;
    for i = 1:size(M_old,1)
        M(i,:) = M(i,:) + NR(i)/size(M_old,2);
    end
end