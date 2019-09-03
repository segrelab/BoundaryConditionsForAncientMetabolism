function [network] = output2network(out,network)
getRxns = @(y) unique(cellfun(@(z) z(1:6), y.input.rxns(~~(y.out.y)),'uni',0));
getMets = @(y) {y.input.compounds(find(y.out.x)).cid};
getRxn_half = @(y) (cellfun(@(z) z, y.input.rxns(~~(y.out.y)),'uni',0));
getRxn_id = @(y) (cellfun(@(z) z(1:6), y.input.rxns(~~(y.out.y)),'uni',0));

% append cid iteration

iter = cellfun(@(x) min(find(x)), num2cell(out.out.X(~~(out.out.x),:),2));
iter_y = cellfun(@(x) min(find(x)), num2cell(out.out.Y(~~(out.out.y),:),2));

%

mets = getMets(out);
rxns = getRxns(out);
rxn_half = getRxn_half(out);
rxn_id = getRxn_id(out);

% issue is here: the reactions intersected are keep half reactions which
% may come online later if the reaction is reversible.
%rx inter makes no snse right now

[~,r,r2] = intersect(network.rxns,rxns);
[~,c,c2] = intersect({network.compounds.cid},mets);

network.S = network.S(c,r);
network.rxns = network.rxns(r);
network.rxn_half = rxn_half;
network.rxn_half_iter = iter_y;
network.compounds = network.compounds(c);
network.met_iter = iter(c2);


rmap = containers.Map(network.rxns,1:length(network.rxns));
[rn,i,j]=unique(rxn_id);
rn_iter = arrayfun(@(k) min(iter_y(j==k)),1:length(rn));
k = cell2mat(cellfun(@(x) rmap(x),rn,'uni',0));
network.rxn_iter = rn_iter(k);


%network.rxn_iter = iter_y(r2);
if isfield(network,'deltaG')
    network.deltaG = network.deltaG(r);
end

if isfield(network,'deltaGerr')
    network.deltaGerr = network.deltaGerr(r);
end

end