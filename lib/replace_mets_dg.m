function [network,coupled_rxns,meta] = replace_mets_dg(network,old,new,ddG)

% find reactions with only couple old and relace with new
if iscell(old)
    map = containers.Map({network.compounds.cid},1:length(network.compounds));
    old = cellfun(@(x) map(x), old,'uni',1);
    new = cellfun(@(x) map(x), new,'uni',1);
end

coupled_rxns = find(sum(~~network.S(old,:)) == 2 & sum(network.S(old,:)) == 0);

for i = 1:length(coupled_rxns)
    s = network.S(old,coupled_rxns(i));
     network.S(new,coupled_rxns(i)) = s;
     network.S(old,coupled_rxns(i)) = zeros(2,1);
     network.deltaG(coupled_rxns(i)) = network.deltaG(coupled_rxns(i)) + -s(1)*ddG;
     meta.s{i} = s;
end

meta.old = old;
meta.new = new;
meta.coupled_rxns = coupled_rxns;
end






