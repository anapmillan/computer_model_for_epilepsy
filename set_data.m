seed = dlmread(seed_filename);
wij = dlmread(network_filename);
load(pattern_filename,'pattern');

nSOZ = numel(seed);
nrois = size(wij,1); %Numer of ROIs

fprintf('Output directory: %s\n', pout)
fprintf('Figures directory: %s\n', pout_figs)

%%
BNV_plots = 0;
extra_plots = 0;

%% Set some values
% Values for beta and gamma (only used if needed according to to w_dyn)
% Can be changed to admit arrays (for testing different values)
beta_val = 0;
gamma_val = 0;
% Max simulation time (only for beta, SIR)
tmax0 = 1e3;

beta = 0;
gamma = 0;
tmax = 0;

%%
n_rho_v = numel(rho_v);
n_beta_v = numel(beta_v);

%% Save all data in variable
create_data_variable;
