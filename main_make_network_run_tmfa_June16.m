
clear
% perform network expansion using defined variables

tic;
%auto = {'C00001','C00011','C00282','C00283','C00288','C00342','C00343','C14818','C14819','C00080','C00138','C00139'};%


%% make a network expansion model using autotrophic network
% H20, CO2, H2S, H2, HCO3, H
auto = {'C00001','C00011','C00282','C00283','C00288','C00080'};
params.CofactorPotential = -220;
params.FeS= true;
params.met_ub = 1e-1;
params.met_lb = 1e-6;
params.keepNanRxns = false;
params.base = 'balanced';
params.CofactorPotentialStoich = 1;
params.ammonia = false;
params.pH = '7.0';
params.Thioester = true;
params.seed = auto;
params.T = 50;
out = run_tne(params);


% build fba model
model = out.irrModel;
model.rxns = model.rxns';
% set the biomass equation
m = addBiomassToProtometabolsim(model,0.2,10);
fparams.met.lb = params.met_lb;
fparams.met.ub = params.met_ub;
fparams.deltag_err_bound = 2;
fparams.T = params.T;

% set the environmental constraints
reductants = 'C00282';
oxidants = 'C00080';
media = {'C00001','C00011',oxidants,reductants};
m2 = setMediaForProtometabolism(m,media',reductants,1);
% run thermodynaimc metabolic flux analysis
o2 = thermoMFA_err(m2,fparams);
t = toc;

