function issuers_int = roundIssuers(issuers, numIssuers)

issuers_int = round(issuers); 
diff = numIssuers-sum(issuers_int); 
if diff ~= 0 
    [~, idx] = sort(abs(issuers-issuers_int), 'descend');
    for i = 1:abs(diff)
        if diff > 0
            issuers_int(idx(i)) = issuers_int(idx(i)) + 1;
        else
            issuers_int(idx(i)) = issuers_int(idx(i)) - 1;
        end
    end
end 
