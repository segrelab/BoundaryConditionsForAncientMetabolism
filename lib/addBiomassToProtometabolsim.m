function [model] = addBiomassToProtometabolsim(model,phi_L,n)


thioester = 'C00024';
thiol = 'C00010';
Rgroup = 'C00033';

reductant = 'C99999';
oxidant = 'C99998';
lipid = 'C00249';



kas = {'C00022';'C00026';'C00036';'C00048';'C00141';'C00168';'C00233';'C00671'};
mm = [88.016 146.0215 132.0054 74.0004 116.05 104.11 130.063 130.063];
[~,idx_ka]=intersect(model.mets,kas);

catalytic_fraction = 1-phi_L;
s_ka = catalytic_fraction/length(kas) * 1e3 ./mm;

[~,idx_lipid]=intersect(model.mets,lipid);
mm = 256.24;
lipid_fraction = phi_L;
s_lipid = catalytic_fraction/length(idx_lipid) * 1e3 ./mm;

% compute electron demand
s_e = 2*sum(s_ka);

% compute thioester demand
s_t = (n-1)/n * sum(s_ka);

mets = [kas',lipid,thioester,thiol,Rgroup,reductant,oxidant,'Biomass'];
%mets = [kas',lipid,thioester,thiol,reductant,oxidant,'Biomass'];

w = [-s_ka,-s_lipid,-s_t,s_t,s_t,-s_e,s_e,1];
%w = [-s_ka,-s_lipid,-s_t,s_t,-s_e,s_e,1];

model = addReaction(model,'biomass_reaction',mets,w,0,0,1000,0);
model = addReaction(model,'biomass_exchange',{'Biomass'},-1,0,0,1000,1);








end