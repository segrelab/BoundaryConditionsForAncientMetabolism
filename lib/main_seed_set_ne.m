function [out] = main_seed_set_ne(input,s)


%[~,~,s]=intersect(seed,{input.compounds.cid});

x0  = zeros(length(input.compounds),1);
x0(s) = 1;
out = netExp(input.R,input.P,x0,input.b);
end
% p = out.x;
% p(9) = 1;
% out_p = netExp(input.R,input.P,p,input.b);
% 
% 
% 
% load('LUCApedia_FeS_ribozyme.mat', 'iron_sulfur')
% rxns = unique(cellfun(@(x) cell2mat(regexp(x,'R\d*(?=\_*)','match')),input.rxns(~~out.y),'uni',0));
% rxns_p = unique(cellfun(@(x) cell2mat(regexp(x,'R\d*(?=\_*)','match')),input.rxns(~~out_p.y),'uni',0));
% rxns_p = setdiff(rxns_p,rxns);