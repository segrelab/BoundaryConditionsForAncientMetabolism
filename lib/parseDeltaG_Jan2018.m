function [dG,dGerr] = parseDeltaG_Jan2018(dGstruct,rns)

g = cellfun(@(x) dGstruct(find(strcmp({dGstruct.rn},x))),rns,'uni',0);
%err = arrayfun(@(x) dGstruct([dGstruct.pg.rn]==x),rns,'uni',0);

for i = 1:length(g)
    if isempty(g{i})
        dG(i) = NaN;
        dGerr(i) = NaN;
    else
        dG(i) = g{i}.deltaG;
        dGerr(i) = g{i}.err;
    end
    
end
end