function [i] = findExchangeReactionsIrrModel(irrModel,met)


    met = ['EX_',met,'_b'];
    i =  find(strcmp(irrModel.rxns,met));



end