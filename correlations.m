function corr = correlations(probDef, S)
    %INPUT
    % probDef: problem default
    % S: size adjustment
    %OUTPUT
    % corr: correlation

    corr = 0.12*(1-exp(-50*probDef))/(1-exp(-50)) + ...
           0.24*(1-(1-exp(-50*probDef))/(1-exp(-50))) - ... 
           0.04*(1-(S-5)/45);  