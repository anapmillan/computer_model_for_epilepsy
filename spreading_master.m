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
%% Model parameters

optimization_options = {'SL','parallel'}; %options for the optmization and dynamics
    % SL (slowlimit): as before, assumes beta--> 0 and optimizes kappa (worst
    %   fit, fastest)
    % parallel: goes over beta and kappa together (best fit, slowest)
w_opt = optimization_options{2}; %select desired option here
 
%% Data files & directories

seed_filename = 'data/seed_example.txt'; % File with seed data (read seed)
network_filename = 'data/network_example.txt'; % File with network data (read wij)
pattern_filename = 'data/pattern_example.mat'; % File with pattern data (load pattern)
% Pattern is a structure indicating the seizure pattern. It has fields:
%   pattern.trois:      (n_sampled_rois x 1) array indicating all sampled ROIs
%   pattern.all_erois:  array indicating the indices of sampled active ROIs
%   pattern.order:      array indicating the activation order of active sampled
%                       ROIs. Size = (n_sampled_active_rois x 1)
%   pattern.ntrois:     total number of sampled ROIs


pout = 'results'; %DO NOT add the final /
pout_figs = 'figures';

w_conn = 'W'; % B for binary, W for weighted

%% Data for loop:
scaling_rho = 'log'; %log, linear
if strcmp(scaling_rho,'log')
    rho_v = logspace(-2,-0.5,10);
else
    rho_v = 0.05:0.05:0.5;
end

scaling_beta = 'log'; 
if strcmp(scaling_beta,'log')
    beta_v = logspace(-5,-0.5,10); 
else
    beta_v = [0.05:0.05:0.5];
end

% Number of iterations (rec: 1e3 - 1e4 (1e2 ok for testing))
nruns = 1e2;

%% Name tag for the gen files
name_tag = sprintf('example_%s_%s_%s',...
                    w_opt, w_conn, nruns);
                               
%% Tags to switch things on and off
 
d_netw = 1; %get netw, threshold, etc
d_matrix_plots = 0; %some nice netw plots (depends on d_netw)

d_dyn = 1; % Run spreading
d_dyn_plots = 0; %Some plots of spreading

d_act_plot = 0; %Better set to 0 for loop (would create too many figures...)
d_corr_i_plot = 0;

d_loop_plot = 1; %Plot final figure
do_verborrea = 1; 
% Control the level of output that you want
% 2 will output several details of the lopps, 1 only the main things, 0
% none

% Some more basic parameters & house keeping
save_spreading_data = 0;
set_data;


%% START RUNNING PROGRAM ACCORDING TO OPTIMIZATION OPTION

if strcmp(w_opt,'slowlimit') 
    model_data.w_dyn = 'SI'; 
    model_data_final = model_data;  
    
    % START LOOPS on: rho (threshold)
    rho_vec_th = zeros(n_rho_v,5);
    fprintf('\nLoop to find best threshold (slow limit)\n')
    
    for i_rho = 1:n_rho_v
        rho = rho_v(i_rho);

        %% Get network
        aij = threshold_matrix_rho(wij,rho,nrois);
        rkm = sum(aij(:)>0)/(nrois*(nrois-1)); %Actual network density

        if (do_verborrea>1 | (do_verborrea & mod(i_rho,5)==0))
            fprintf('   Network density: %.2f (%d/%d)\n', rho, i_rho, n_rho_v)       
        end
        %% Run spreading, measure activation, compare with PET
 
        model_data.aij = aij;
        [~, rho_vec_th(i_rho,:)] = ...
            run_dyn(model_data, pattern); 
      
    end
    
    %% Find best threshold  
    [cmax,jmax] = max(rho_vec_th(:,5));
    rho_max = rho_v(jmax);    
    fprintf('Best correlation c= %.3f found for th = %.4f\n', ...
        cmax, rho_max)
    
    %Correlation figure:
    n_fig_SL = sprintf(...
        '%s/SL_rho_fit_%s.png', pout_figs, nruns);
    make_SL_figure; %figure tag: f_SL
    print(f_SL,'-dpng',n_fig_SL);
    
    % Prepare for the final plots:
    final_rho = rho_max;    
    tsave = struct();
    tsave.rho_vec_th = rho_vec_th;
    tsave.model_data = model_data;
    tsave.cmax_th = cmax; 
elseif strcmp(w_opt,'parallel')
    fprintf('\nParallel mode: fitting threshold and spreading rate together\n')
    
    model_data.tmax = 1e3;
    model_data.w_dyn = 'beta';
    model_data_final = model_data;
    rho_vec_parallel = zeros(n_rho_v, n_beta_v,5);
    
    tic
    %Loop on thresholds
    for i_rho = 1:n_rho_v
        
        %Threshold matrix to new density
        rho = rho_v(i_rho);
        aij = threshold_matrix_rho(wij,rho,nrois);
        rkm = sum(aij(:)>0)/(nrois*(nrois-1)); %Actual netw density
        model_data.aij = aij;          
        if do_verborrea
            fprintf('\nNetwork density: %.2f (%.2f) %d/%d\n', rho, rkm, i_rho, n_rho_v)
        end
        
        %Loop on spreading rates
        for i_beta = 1:n_beta_v
            if (do_verborrea>1 | (do_verborrea & mod(i_beta,5)==0))
                fprintf('    Spreading rate: %.3f\n', beta_v(i_beta))
            end
            model_data.beta = beta_v(i_beta);
            [~, rho_vec_parallel(i_rho,i_beta,:)] = run_dyn(model_data, pattern); 
        end
        toc
    end
    %% Find best parameters
    aux = squeeze(rho_vec_parallel(:,:,5));
    [cmax_parallel, idx] = max(aux(:));
    [imax, jmax] = ind2sub(size(aux), idx);
    fprintf(...
        'Best fit c = %.4f obtained for rho = %.3f, beta = %.3f\n',...
        cmax_parallel, rho_v(imax), beta_v(jmax))
    
    %Make plots 
    n_fig_parallel = sprintf(...
        '%s/parallel_fit_%d.png', pout_figs, nruns);    
    make_parallel_figure; %f_parallel
    print(f_parallel, '-dpng', n_fig_beta);
    
    %Prepare for the final plots
    model_data_final.beta = beta_v(jmax);
    rho_max = rho_v(imax);
       

    tsave = struct();
    tsave.rho_vec_parallel = rho_vec_parallel;
    tsave.cmax_parallel = cmax_parallel;
         
else
    error('The optimization option makes no sense!\n Bye...\n')
end


tsave.model_data_final = model_data_final; 
tsave.model_data = model_data; 
save(sprintf('%s/master_savefile_optimization_%s.mat', model_data.pout_figs, ...
    model_data.name_tag),'tsave')


%% Plot details of best fit
model_data_final.d_dyn_plots = 1;
model_data_final.d_act_plot = 1;
model_data_final.d_corr_i_plot = 1;
model_data_final.name_tag = [model_data_final.name_tag, '_best_fit'];

aij = threshold_matrix_rho(wij,rho_max,nrois);
model_data_final.aij = aij;
[tot_pob, rho_final] = ...
    run_dyn(model_data_final, pattern);

