function [acu_tot_pob, corr_vec] = run_dyn(model_data, pattern)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If you use this code, then please cite:
%  1.- Ana P Millan et a., "Epidemic models characterize seizure propagation 
%      and the effects of epilepsy surgery in individualized brain networks 
%      based on MEG and invasive EEG recordings." medRxiv (2021).
%  2.- Ida Nissen et al. "Optimization of epilepsy surgery through virtual 
%      resections on individual structural brain networks." 
%      Scientific Reports 11.1 (2021): 1-18.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [acu_tot_pov, rho_vec, pval_vec] = run_dyn(model_data, pattern)
% runs spreading_func to obtain spreading, calculates acu_tot_pob and compares
% with pattern
% INPUTs:
%   model_data: structure variable with all data for the model simulation
%               generated according to create_data_variable
%   pattern:    structure variable with pattern information
%               


%
% OUTPUTs:
%   acu_tot_pob:    cumulative probability distribution of infection orders
%                   (n_rois x n_orders) matrix                    
%   corr_vec:       correlation data
%                   array including: [
%                        pcp: pearson's correlation
%                        coefficient, frac_act: fraction of sampled rois active in
%                           model and data; 
%                        frac_eq: fraction of sampled rois in the
%                        same state (active/inactive) in model and data;
%                        n_act_model: fraction of active rois in the model;
%                        correlation metrix: pcp*frac_eq]

    %% Run spreading mdoel
    tot_pob = spreading_func(model_data.aij, model_data.seed, model_data.nruns, ...
        model_data.name_tag, model_data.pout, model_data.w_conn,...
        model_data.w_dyn, model_data.beta, model_data.tmax, ...
        model_data.save_data, model_data.d_dyn_plots); 
    % close all;

    acu_tot_pob = zeros(size(tot_pob));
    for i = 1:size(tot_pob,2)
        acu_tot_pob(:,i) = sum(tot_pob(:,1:i),2);
    end
    if model_data.d_act_plot; make_d_act_plot; end


    %% Compare with pattern

    aux = corr_model_data(acu_tot_pob, pattern);
    aux = cell2mat(aux);
    corr_vec = zeros(1,5);
    corr_vec(1:4) = aux;
    corr_vec(5) = aux(1)*aux(3);
