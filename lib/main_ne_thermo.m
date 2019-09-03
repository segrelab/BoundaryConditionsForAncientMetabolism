function [output] = main_ne_thermo(input,seed,ub,lb,keep_nan,temp)


RT = 0.008309424 * (273.15+temp);
% remove h2o and hyrodgen ions from this calculation
[~,mm]=intersect({input.compounds.cid},{'C00001','C00080'});

Sreact = abs(input.S.*input.R);
Sprod = abs(input.S.*input.P);

Sreact(mm,:) = [];
Sprod(mm,:) = [];

k = RT*log((lb.^sum(Sprod))./(ub.^sum(Sreact))) + input.deltaG > 0;

if ~keep_nan
   k = k | isnan(input.deltaG);
end
        


    input.P(:,k) = [];
    input.R(:,k) = [];
    input.rxns(k) = [];
    input.matched(k) = [];
    input.deltaG(k) =[];
    input.deltaGerr(k) = [];
    input.b(k) = [];
    
    % remove 
    c = ~sum(input.R') & ~sum(input.P');
    input.compounds(c) = [];
    input.R(c,:) = [];
    input.P(c,:) = [];






if isempty(seed)
    out = main_seed_ne(input);
    
else
    [~,seed] = intersect({input.compounds.cid},seed);
    out = main_seed_set_ne(input,seed);
end

output.out = out;
output.input = input;
end

