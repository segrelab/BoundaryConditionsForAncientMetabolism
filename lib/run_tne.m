function [out,network] = run_tne(params) 

CofactorPotential = params.CofactorPotential;
pH = params.pH;
Temp = params.T;
seed = params.seed;
FeS = params.FeS;
Thioester = params.Thioester;
ub = params.met_ub;
lb = params.met_lb;
keepNanRxns = params.keepNanRxns;

switch params.base
    case 'consistent'
       load('networks/network_consistent','network')
    case 'balanced'
        load('network_filtered.mat', 'network');
    case 'no_mutlistep'
        load('networks/network_consistent_noMultiStep_Feb26','network');
    case 'full'
        load('networks/KEGG_network_full.mat');
    case 'p_balanced'
        load('networks/KEGG_network_P_balanced.mat')
    case 'balanced_no_multistep'
        %load('networks/network_balanced_noMultiStep_Feb26','network');
        load('networks/network_balanced_noMultiStep_March3','network');
end


%base_network = params.base_network;
%load('network_filtered.mat', 'network');
%load('networks/network_consistent','network')


og = network;

if ~isnan(CofactorPotential)
    % remove flavins and nicotinamide transhydrogenases
    [~,ii]=intersect(network.rxns,{'R10159';'R01195';'R00112';'R09520';'R09748';'R05705';'R05706';'R09662';'R09750'});
    network = removeReactionsFromNetwork(network,network.rxns(ii));
end


%load('networks/KEGG_network_full.mat');
%load('networks/KEGG_network_P_balanced.mat');
%[network,usesCoa] = subs_pantetheine_for_coA(network);
%seed = importdata('seeds/cononical_seed.csv');

%seeds.cho = {'C00001','C00011','C00288','C00033','C00058'};
%seeds.s = {'C00283'};
%seeds.n = {'C00697','C00014'};

% delete reaction that degrades pantetheine
%network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R02973')));
network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R10747')));
%network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R02972')));
% for glutathione usage
network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R00494')));
network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R01262')));
network = removeReactionsFromNetwork(network,find(strcmp(network.rxns,'R03916')));

% find benzoate deg. reactions
% data = urlread('http://rest.kegg.jp/link/rn/map00362');
% benzoate = regexp(data,'(?<=rn:)R\d*','match');
% [~,benzoate]=intersect(network.rxns,benzoate);
% network = removeReactionsFromNetwork(network,benzoate);

% find all reactions that a uniquely part of degradative pathways
% load('/Users/joshuagoldford/Dropbox/phd/data/KEGG/modules/modules.mat', 'mods');
% deg_mods = find(not(cellfun(@isempty,regexp(mods.name,'degradation'))));
% nd_mods = find((cellfun(@isempty,regexp(mods.name,'degradation'))));
% torem =setdiff(find(sum(mods.M(deg_mods,:)) > 0),find(sum(mods.M(nd_mods,:)) > 0));
% [~,torem]= intersect(network.rxns,mods.rxns(torem));
% network = removeReactionsFromNetwork(network,torem);
% disp('removed degradation modules...');

% remove iron sulfur reactions
if ~FeS
    %load('LUCA_RXNS.mat');
    iron_s_rxns = load('FeS_rxns_EC_KEGG');
    %[~,fes_r] = intersect(network.rxns,data(4).reactions);
    [~,fes_r] = intersect(network.rxns,iron_s_rxns.FeS);
    network = removeReactionsFromNetwork(network,fes_r);
    disp('removed iron sulfer reaction...');
end

%deltaG = load('deltaG_kegg');
addpath reaction_free_energy
deltaG = importKeggDeltaG(['reaction_free_energy/kegg_reactions_CC_ph',pH,'.csv']);
[dG,dGerr] = parseDeltaG_Jan2018(deltaG,network.rxns);
network.deltaG = dG';
network.deltaGerr = dGerr';
disp('added free energy estimates to the network...');


phos.mono = 'C00009';
phos.di = 'C00013';
phos.tri = 'C00536';
phos.ddG = 15.5;


X = zeros(10,length(network.rxns))';
% repace all NTP - NDP coupled reactions with P-P_i --> P_i
% make a P-P + P pair for ATP, ADP
[np,c] = replace_mets_dg(network,{'C00002','C00008'},{phos.di,phos.mono},phos.ddG);
X(c,1) = 1;
% replace GTP, GDP
[np,c] = replace_mets_dg(np,{'C00044','C00035'},{phos.di,phos.mono},phos.ddG);
X(c,2) = 1;
% replace CTP, CDP
[np,c]  = replace_mets_dg(np,{'C00063','C00112'},{phos.di,phos.mono},phos.ddG);
X(c,3) = 1;

% replace UTP, UDP
[np,c]  = replace_mets_dg(np,{'C00075','C00015'},{phos.di,phos.mono},phos.ddG);
X(c,4) = 1;

% replace ITP, IDP
[np,c]  = replace_mets_dg(np,{'C00081','C00104'},{phos.di,phos.mono},phos.ddG);
X(c,5) = 1;

% for diphosphate transfers
% replace ATP, AMP
[np,k]  = replace_mets_dg(np,{'C00002','C00020'},{phos.tri,phos.mono},0);
X(k,6) = 1;

% replace GTP, GMP
[np,k]  = replace_mets_dg(np,{'C00044','C00144'},{phos.tri,phos.mono},0);
X(k,7) = 1;

% replace CTP, CMP
[np,k]  = replace_mets_dg(np,{'C00063','C00055'},{phos.tri,phos.mono},0);
X(k,8) = 1;

% replace UTP, UMP
[np,k]  = replace_mets_dg(np,{'C00075','C00105'},{phos.tri,phos.mono},0);
X(k,9) = 1;

% replace ITP, IMP
[np,k]  = replace_mets_dg(np,{'C00081','C00130'},{phos.tri,phos.mono},0);
X(k,10) = 1;

network = np;

disp('subs phosphate coupling for primitive coupling...');


%% remove cofactor coupling
%CofactorPotential = -240;
%[network,redox] = replaceRedoxCouplings(network,CofactorPotential);
%[network,redox] = replaceCoenzymeCouple(network,'C00003','C00004',CofactorPotential);

%% old replacement of coenzymes
% if ~isnan(CofactorPotential)
%     [network,redox] = replaceRedoxCouplings(network,CofactorPotential);
%     seed = union(seed,{'C99998';'C99999'});
% end
    

%% new replacement of coenzymes
if ~isnan(CofactorPotential)
    % replace NAD/FAD/Frdxn/Trdxn/CoQ
    load('half_potentials','e0')
    ph_num = str2num(pH);
    CofactorPotentialStoich = params.CofactorPotentialStoich;
    % swap nad rxns
    n0 = network;
    potential = e0.nad([e0.nad.pH] == ph_num).E0;
    [network,rrn] = replaceCoenzymeCouple(network,'C00003','C00004',[1 1],CofactorPotentialStoich,potential,CofactorPotential);
    disp(['removed NAD redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

    % swap nadp rxns
    potential = e0.nadp([e0.nadp.pH] == ph_num).E0;
    network = replaceCoenzymeCouple(network,'C00006','C00005',[1 1],CofactorPotentialStoich,potential,CofactorPotential);
    disp(['removed NADP redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

     % swap fad rxns
    potential = e0.fad([e0.fad.pH] == ph_num).E0;
    network = replaceCoenzymeCouple(network,'C00016','C01352',[1 1],CofactorPotentialStoich,potential,CofactorPotential);
    disp(['removed FAD redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

    % swap coq rxns
    %potential = e0.coq([e0.coq.pH] == ph_num).E0;
    %network = replaceCoenzymeCouple(network,'C00390','C00399',[1 1],potential,CofactorPotential);
    %disp(['removed CoenzymeQ redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

    % swap glutathione rxns
    %potential = e0.glt([e0.glt.pH] == ph_num).E0;
    %network = replaceCoenzymeCouple(network,'C00127','C00051',[1 2],potential,CofactorPotential);
    %disp(['removed Glutathione redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

    % swap thioredoxin rxns
    %potential = e0.txn([e0.txn.pH] == ph_num).E0;
    %network = replaceCoenzymeCouple(network,'C00343','C00342',[1 1],potential,CofactorPotential);
    %disp(['removed thioredoxin redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

    % swap ferrodoxin rxns
    %potential = e0.fxn([e0.fxn.pH] == ph_num).E0;
    %network = replaceCoenzymeCouple(network,'C00139','C00138',[1 1],potential,CofactorPotential);
    %disp(['removed ferrodoxin redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);
    seed = union(seed,{'C99998';'C99999'});
end




%rns.redox = [network.rxns(redox.nad);network.rxns(redox.nadp);network.rxns(redox.fad)];
%disp(['removed redox coupling and replacing with a generic coenzyme with a potential of ', num2str(CofactorPotential), ' mV']);

% convert to input structure
network = removeUnusedMetabolites(network);
[network,irrev] = findIrrevRxns(network);
input = network2inputStructure_dgerr(network,irrev);
input = append_input_reversible_matched(input);
disp('constructed input structure...');

%input = add_free_energy_to_input(input);
%dg_cutoff = 0:1:60;
%pan = 'C00831';
%acp = 'C00227';
disp('...running network expansion...');


%out = main_ne_thermo(input,[seeds.cho';'C99998';'C99999';pan],1e-2,1e-6,true);
%out = main_ne_thermo(input,[seeds.cho';seeds.s';'C99998';'C99999';pan],1e-1,1e-6,false);


if Thioester
    seed = union(seed,'C00010');
    seed = union(seed,'C00229');
    seed = union(seed,'C00051');
end
%seed = union(seed,{'C99998';'C99999'});

ne_res = main_ne_thermo(input,seed,ub,lb,keepNanRxns,Temp);
disp('...Finished!');

% constructing model from net exp
[model,irr] = constructModelFromNetworkExpansion2(ne_res,network,seed);
%getRxns = @(out) unique(cellfun(@(x) x(1:6),out.input.rxns(~~out.out.y),'uni',0));
%rxns = getRxns(out);
n = output2network(ne_res,network);


% remove highly connected metabolites
[Adj,k,R] = Stoic2Adj(n.S,n.rxn_iter,0.99);
Adj = R;

cids = {n.compounds.cid};
compounds = n.compounds;
compounds(k) = [];

met_iter = n.met_iter;
met_iter(k) = [];


A = cell2mat(arrayfun(@(x) abs(x-met_iter)',met_iter,'uni',0));
adj = Adj.*double(A==1);
% save nodes for generic d3 usage
id = {compounds.cid};
iter = num2cell(met_iter)';
name = reducedNameSet({compounds.name});
compounds(end-1).elements = zeros(118,1);
compounds(end).elements = zeros(118,1);

numC = cellfun(@(x) x(6), {compounds.elements},'uni',0);
%iscoa = cellfun(@(z) ~isempty(z),regexp({compounds.name},'CoA'),'uni',0);
iscoa = cellfun(@(z) ~isempty(z),regexp({compounds.name},'CoA|\[acp\]'),'uni',0);
iscoao = cellfun(@(z) ~isempty(z),regexp({compounds.name},'CoA'),'uni',0);

for i = 1:length(iscoao)
    if iscoao{i}
        numC{i} = numC{i} -21;
    end
end

nodes = struct('id',id,'iter',iter,'cpdname',name,'numCarbons',numC,'CoA',iscoa);
nodes = struct2table(nodes);
E = adj2edgeL(adj);
[~,e]=unique(cell2mat(cellfun(@sort,num2cell(E(:,1:2),2),'uni',0)),'rows');
E = E(e,:);

edge = table();
edge.source = E(:,1);
edge.target = E(:,2);
edge.value = E(:,3);
graph.nodes = nodes;
graph.edges = edge;


model.mets = n.compounds;

out.revModel = model;
out.irrModel = irr;
out.graph = graph;
out.params = params;
out.network = n;

% decrement index for json
%E(:,1:2) = E(:,1:2)-1;
%links = struct('source',cellfun(@(x) x,num2cell(E(:,1),2),'uni',0),'target',cellfun(@(x) x,num2cell(E(:,2),2),'uni',0),'value',num2cell(E(:,3),2));
%links = struct('source',cellfun(@(x) nodes(x+1).id,num2cell(E(:,1),2),'uni',0),'target',cellfun(@(x) nodes(x+1).id,num2cell(E(:,2),2),'uni',0),'value',num2cell(E(:,3),2));
% 
% graph.nodes = nodes';
% graph.links = links;
% %savejson('graph',graph,'ne_CHOThioester_760mvThermo_Jan13.json');
% 
end
