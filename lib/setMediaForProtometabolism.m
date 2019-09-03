function [model] = setMediaForProtometabolism(model,media_set,limiting_metabolite,limitingFlux)



g = cellfun(@(g) findExchangeReactionsIrrModel(model,g),model.mets,'uni',0);
model.ub(cell2mat(g)) = 0;

uub = 1e3;
media_set = intersect(model.mets,media_set);
z = cellfun(@(x) findExchangeReactionsIrrModel(model,x),media_set);
model.ub(z) = uub;
model.ub(findExchangeReactionsIrrModel(model,limiting_metabolite)) = limitingFlux;


end