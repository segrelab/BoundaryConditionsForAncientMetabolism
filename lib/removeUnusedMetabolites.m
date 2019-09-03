function [network] = removeUnusedMetabolites(network)
%% Author: J. Goldford
% 

S = ~~network.S;

to_remove = sum(S,2) == 0;
network.S(to_remove,:) = [];
if isfield(network,'compounds')
network.compounds(to_remove) = [];
else
    network.mets(to_remove) = [];
    network.b(to_remove) = [];
end

end





