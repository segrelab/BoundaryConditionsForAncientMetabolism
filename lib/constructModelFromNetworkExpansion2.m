function [model,irr] = constructModelFromNetworkExpansion2(out,network,seed)


%load('network_filtered.mat', 'network');
%seed = importdata('seeds/cononical_seed.csv');
%network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R10747')));
%coa = 'C00010';
%hetero = main_ne_network_dg(network,[seed;'C00010'],[],[]);

n0 = network;
hetero = out;

%% can we synthesize all components in network using delta G cutoff

getRxns = @(y) unique(cellfun(@(z) z(1:6), y.input.rxns(~~(y.out.y)),'uni',0));
getMets = @(y) {y.input.compounds(find(y.out.x)).cid};

map.c = containers.Map({n0.compounds.cid},1:length(n0.compounds));
map.r = containers.Map(n0.rxns,1:length(n0.rxns));
map.cids = @(y) cellfun(@(x) map.c(x),y);
map.rxns = @(y) cellfun(@(x) map.r(x),y);


model.rxns = getRxns(hetero);
model.mets = getMets(hetero);
model.S = n0.S(map.cids(model.mets),map.rxns(model.rxns));
model.rev = ones(length(model.rxns),1);
model.lb = -1e4.*ones(length(model.rxns),1);
model.ub = 1e4.*ones(length(model.rxns),1);
model.c = zeros(length(model.rxns),1);
model.b = zeros(length(model.mets),1);
model.deltaG = network.deltaG(map.rxns(model.rxns));
model.deltaGerr = network.deltaGerr(map.rxns(model.rxns));

seed = intersect(seed,model.mets);
k = cellfun(@(x) find(strcmp(x,model.mets)),seed);
lb = -10;
ub = 1e4;
for i = 1:length(model.mets)
    %if ismember(i,k)
        model = addExchangeRxn(model,model.mets(i),lb,ub);
    %else
       % model = addExchangeRxn(model,model.mets(i),0,ub);
    %end 
end
model.deltaG = [model.deltaG',nan(1,length(model.rxns)-length(model.deltaG))]';
model.deltaGerr = [model.deltaGerr',nan(1,length(model.rxns)-length(model.deltaGerr))]';

reverse = find(model.lb < 0 & model.ub > 0);



irr.S = [];
irr.rxns = {};
irr.mets = model.mets;
irr.lb = [];
irr.ub = [];
irr.c = [];
irr.b = model.b;
irr.deltaG = [];
irr.deltaGerr= [];
z = 1;
for i= 1:length(model.rxns)
    if ismember(i,reverse)
        irr.S = [irr.S model.S(:,i)];
        irr.S = [irr.S -model.S(:,i)];
        irr.rxns = [irr.rxns,[model.rxns{i},'_f']];
        irr.rxns = [irr.rxns,[model.rxns{i},'_b']];
        irr.lb = [irr.lb 0 0];
        irr.ub = [irr.ub model.ub(i) -model.lb(i)];
        irr.deltaG = [irr.deltaG model.deltaG(i) -model.deltaG(i)];
        irr.deltaGerr = [irr.deltaGerr model.deltaGerr(i) model.deltaGerr(i)];

    else
        irr.S = [irr.S model.S(:,i)];
        irr.rxns = [irr.rxns,[model.rxns{i},'_f']];
        irr.lb = [irr.lb 0];
        irr.ub = [irr.ub model.ub(i)];
        irr.deltaG = [irr.deltaG model.deltaG(i)];
        irr.deltaGerr = [irr.deltaGerr model.deltaGerr(i)];

    end
end
irr.c = zeros(length(irr.rxns),1);
irr.lb = irr.lb';
irr.ub = irr.ub';
irr.mets = irr.mets';
irr.deltaG = irr.deltaG';
irr.deltaGerr = irr.deltaGerr';



%irr.model = convertToIrreversible(model);
% 
% 
% map.irr = containers.Map(hetero.input.rxns,1:length(hetero.input.rxns));
% map.irr_rxns = @(y) cellfun(@(x) map.irr(x),y);
% 
% irr.model.delta_g = hetero.input.deltaG(map.irr_rxns(irr.model.rxns));
% 
% 
% irr.model = addExchangeRxn(irr.model,seed,-1000.*ones(length(seed),1),1000.*ones(length(seed),1));
% r = setdiff(irr.model.mets,seed);
% 
% irr.model = addExchangeRxn(irr.model,r,1.*ones(length(r),1),1000.*ones(length(r),1));
% irr.model.b = zeros(length(irr.model.mets),1);

end