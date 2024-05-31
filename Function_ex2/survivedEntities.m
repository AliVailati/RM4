function [survived] = survivedEntities (InitialDistribution, years, M_avg1y)

survived = InitialDistribution;
for ii = 1:years
    survived = survived .* M_avg1y;
end
survived = diag(survived); 
numPopulation = round(sum(survived)); 
survived = roundIssuers(survived, numPopulation);