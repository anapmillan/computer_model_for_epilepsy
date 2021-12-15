function aij = network_func(w_netw, rho, seed, name_tag, d_matrix_plots, pout)
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
% function aij = network_func(netw_file, rho, seed, name_tag, d_matrix_plots, pout)
% netw_file: inputfile including path
% rho: netw density (-1 for MST)
% seed: seed rois (for plots)
% name_tag: tag string for file names
% d_matrix_plots: switch to control output (adj matrix plots)
% pout: directory for outputs (if plots on)

    if strcmp(w_netw,'euclidean') & strfind(name_tag,'W')
        fprintf('You asked for weighted eucliean network must be binary\n')
        fprintf('So Im creating a binary network\n')
        fprintf('For a weighted Euclidean network a rule must be defined first: \n')
        fprintf('How do the weights scale with the distance?\n')
    end

    %Read or create network
    aux = 0;
    if strcmp(w_netw,'AEC')
        %File storing network (including path)
        netw_file = '/data/analyses/Ana_Deborah/data/AEC_broadband_SCDneg.txt';
    elseif strcmp(w_netw,'AECC')
        netw_file = '/data/analyses/Ana_Deborah/data/export_matrix_SCDneg_AECc_beta.txt';
    else
        netw_file = '../data/SCDneg_distance_matrix.txt';  
    end

    fprintf('Reading from file: %s\n', netw_file)
    if ~isfile(netw_file)
        error('ERROR: Network file does not exits!/n');       
    else %Read network from file
        wij = dlmread(netw_file);
    end
    if strcmp(w_netw, 'exp') | strcmp(w_netw,'euclidean')  
        %We make a trick for the euclidean network: weight the links
        %accoring to exp rule which is based in distance, then threshold to
        %keep stronger ones and binarize (and therefore remove exponential
        %rule)
        %For the binary version it does not matter which scaling we use for
        %the weights, it is just a tool to order them (i.e. it could just be
        % -distance but I'm not sure what my existing functions would do with negative weights only)
        
        
        %Create network with exponential rule
        fprintf('Creating network with exponential rule\n')
        wij = exp(-0.188*wij);
        wij = wij - diag(diag(wij));
        fprintf('Hi, there might have been an error before so that the diagonal was not 0 for the exponential distance network\n')
        fprintf(' I am not sure if it is adressed later on, maybe best to repeat that analysis.\n')
    end
    
%     figure
%     imagesc(wij)
%     colorbar
%     pause
    
    %Make sure it's symmetric and normalize
    wij = 0.5*(wij+wij');
    wij = wij/max(wij(:));

  %  imagesc(wij)
  %  rho
    kmo = rho*size(wij,1);
 %   pause
    if rho>0 
        bij = get_BINadj_B(wij,kmo);
        bij = (bij + bij')>0;
    else
        fprintf('Getting MST\n');
        bij = kruskal_algorithm(wij);
    end
    
    if strfind(name_tag,'W') & ~strcmp(w_netw,'euclidean')
        aij = wij .* bij;
    else 
        aij = bij;
    end
    %aij
    clear wij 
   
    if d_matrix_plots
        rkm = sum(aij(:));
        BNV_plots = 0;
        network_plots_func(aij,name_tag,pout,seed,rkm,BNV_plots);
    end
end