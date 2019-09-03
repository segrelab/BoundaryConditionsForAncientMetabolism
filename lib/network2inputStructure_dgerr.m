function [out] = network2inputStructure_dgerr(network,irrev)
% author: Joshua Goldford
% date: 11-6-2015
%% network2inputStructure produces the input structure from network structure and irreversibility vector
% inputs: network := network structure with fields: mxn'S', mx1'compounds',
%                   nx1 'rxns'
%         irrev   := nx1 vector, 1 if reaction is irreversible 0 otherwise



rxns = network.rxns;

% convert stoichiometrix matrix into nx1 cell
S = num2cell(network.S,1);

% find reversible reactions
rev = find(~irrev);

ir = find(irrev);
deltaG = num2cell(network.deltaG);
deltaGerr = num2cell(network.deltaGerr);

% for each reversible reaction
for i =1:length(rev)
    % convert single stoichiometric vector into both forwar and backward
    % reactions
    S{rev(i)} = [ S{rev(i)},  -S{rev(i)}];
    % update reaction names with _f for forward and _b for backward
    rxns{rev(i)} = {strcat(rxns{rev(i)},'_f'),strcat(rxns{rev(i)},'_b')};
    deltaG{rev(i)} = [deltaG{rev(i)},  -deltaG{rev(i)}];
    deltaGerr{rev(i)} = [deltaGerr{rev(i)},deltaGerr{rev(i)}];

end

% parse S-cell into S-matrix
S = cell2mat(S);
% parse reactions
rxns = [rxns{:}]';

% reactant matrix is 1 is s(i,j) < 0, 0 otherwise
out.S = S;
out.R = double(S < 0);
% product matrix is 1 is s(i,j) > 0, 0 otherwise
out.P = double(S > 0);
% b(j) is the sum of all reatants for reation j
out.b = sum(out.R)';
% append the compound structure from the original network
out.compounds = network.compounds;
% append updated reaction network
out.rxns = rxns;
out.deltaG = cell2mat(deltaG');
out.deltaGerr = cell2mat(deltaGerr');

end


