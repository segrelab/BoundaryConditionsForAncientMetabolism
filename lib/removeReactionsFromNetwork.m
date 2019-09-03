function [network] = removeReactionsFromNetwork(network,rxns)

if ~iscell(rxns);
    network.S(:,rxns) = [];
    network.rxns(rxns) = [];
    
    
    network = removeUnusedMetabolites(network);
else
    [~,rxns] = intersect(network.rxns,rxns);
    network.S(:,rxns) = [];
    network.rxns(rxns) = [];
    network = removeUnusedMetabolites(network);
end
end