function [network,irrev,r] = findIrrevRxns(network)


o2  = find(strcmp('C00007',{network.compounds.cid}));

%make all reactions that produce 02 consume o2
r = find(network.S(o2,:)' > 0);
network.S(:,r) = -network.S(:,r);
irrev = (network.S(o2,:)' < 0);
end