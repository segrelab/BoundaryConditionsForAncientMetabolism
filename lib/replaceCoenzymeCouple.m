function [network,redox_reactions] = replaceCoenzymeCouple(network,oxidant,reductant,stoich,stoichMulti,origPotential,CofactorPotential)
% function replaces nad/nadp and fad redox coenzymes with fake coenzyme
% with fixed standard potential

% input is a network structure and cofactor potential in mV rel. to SHE
ox_generic = 'Oxidant';
red_generic = 'Reductant';
oxidant_str = struct('cid','C99998','formula','R','smiles','R','name','Generic oxidant','num_rgroup',1,'elements',[],'electrons',nan);
reductant_str = struct('cid','C99999','formula','R','smiles','R','name','Generic reductant','num_rgroup',1,'elements',[],'electrons',nan);

% 
ox = find(strcmp({network.compounds.cid},'C99998'));
red = find(strcmp({network.compounds.cid},'C99999'));
S = network.S;
cpds = network.compounds;

if isempty(ox)
    cpds = [network.compounds,oxidant_str,reductant_str];
    ox = length(cpds)-1;
    red = length(cpds);
    % now add to network generic coenzymes
    S = [S;zeros(2,length(network.rxns))];
end


n = 2;
F = 96.485;

ddG = 0.001*n*F*(CofactorPotential - origPotential);


% find all the reactions

couple = {oxidant,reductant};

findMet = @(x) find(strcmp(x,{network.compounds.cid}));
met_idx = cellfun(@(y) findMet(y), couple);


%intersect({network.compounds.cid},nad);
%[~,nadp_idx] = intersect({network.compounds.cid},nadp);
%[~,fad_idx] = intersect({network.compounds.cid},fad);

% make S matrix;
redox_reactions = find(abs(sum(S(met_idx,:))) == diff(stoich) & sum(~~S(met_idx,:)) == sum(stoich));

% for each redox set, determine if coenzyme is the oxidant
[~,oxidant_idx] = intersect({network.compounds.cid},oxidant);

% for each redox set, determine if coenzyme is the oxidant
[~,reductant_idx] = intersect({network.compounds.cid},reductant);



% reduction potentials are defined like OX + 2e- = Red;
% if reaction is oxidant then ddG = 0.001*nF(CofactorPotentia - NAD_p)
% if reactant is reductant, then ddG = -0.001*nF(CofactorPotentia - NAD_p)
redox_dir = S(oxidant_idx,redox_reactions) < 0;
network.deltaG(redox_reactions) = network.deltaG(redox_reactions) + ddG.*(-1).^redox_dir';



% swap coenzyme stoichiometries
S([ox;red],redox_reactions) = stoichMulti.*S(met_idx,redox_reactions);

% remove redox couples from reactions

S(met_idx,redox_reactions) = 0;

% remove redox couples from reactions
network.S = S;
network.compounds = cpds;
end