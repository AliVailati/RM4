function M = removeNR(M_old)
    %INPUT
    % M_old: matrix with NR
    %OUTPUT
    % M: matrix without NR

    %Before removing NR, we have to evaluate if all the row of M_old sum to 100
    
    %Check if the sum of the rows is equal to 100
    if sum(sum(M_old, 2) == 100) ~= size(M_old, 1)
        for ii = 1:size(M_old, 1)
            if sum(M_old(ii, :)) ~= 100
                diff = 100 - sum(M_old(ii, :));
                M_old(ii, end) = M_old(ii, end) + diff;
            end
        end
    end
    %NR is the last column of M_old
    NR = M_old(:,end);
    %Remove NR from M_old
    M_old = M_old(:,1:end-1);
    %Distribute NR on the rows with respect to the weight of each element in the row
    M = M_old;
    for i = 1:size(M, 1)
        weights = M(i, :)./sum(M(i, :));
        M(i, :) = M(i, :) + NR(i)*weights;
    end 
