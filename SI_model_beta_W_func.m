function [order_in,time_in] = SI_model_beta_W_func(aij,seed,beta,nruns,tmax)
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
% function order_in = SI_model_beta_W_func(aij, seed, beta, nruns, tmax);
% SI_model_B_func runs the SI model for a given beta.
% INPUTs:
%   aij:    weighted adjacency matrix
%   seed:   seed regions (array)
%   beta:   spreading rate
%   nruns:  number of repetitions
%   tmax:   maximum number of simulation steps
%
% OUTPUTs:
%   order_in:   probability map of each ROI to become infected at each  order 
%               (not normalized).
%   time_in:    probability map of each ROI to become infected at each step 
%               (not normalized).

    nr_ROIs = size(aij,1);
    frate = 1;  % Fraction of connected ROIs to be infected before dynamics stops
    
    order_in = zeros(nr_ROIs+10,nr_ROIs);
    time_in = zeros(tmax,nr_ROIs);


    % Create Degree array and neighbour matrix
    bij = aij>0;                        % Binary matrix
    ki = sum(bij);
    kmax = max(ki);
    v_neigh = zeros(nr_ROIs,kmax);      % Neighbour matrix (from binary network)
    v_neigh_w = zeros(nr_ROIs,kmax);    % Weighted neihbour matrix (for probabilities)
    for i = 1:nr_ROIs
        v_neigh(i,1:ki(i)) = find(bij(i,:));
        v_neigh_w(i,1:ki(i)) = aij(i,v_neigh(i,1:ki(i)));
    end
 
    % Calculate number of connected ROIs (to stop dynamics when all are
    % infected)
    % First: distance to the seed d_soz
    gij = graph(aij);               
    dij = distances(gij);               % Distance matrix
    if numel(seed)==1
        d_soz = dij(seed,:);
    else        
        d_soz = min(dij(seed,:));
    end
    % Connected ROIs: finite distance to the seed
    con_rois = find(isfinite(d_soz));
    n_con_rois = numel(con_rois);
    
    % Loop over iterations
    parfor ir = 1:nruns        
        v_neigh_in = v_neigh;
        v_neigh_w_in = v_neigh_w;
        
        v_order = zeros(nr_ROIs,1); %Matrix of infection order
        v_time = zeros(nr_ROIs,1);  %Matrix of infection times
        
        % Seed
        vinf = seed';                % Array of infected ROIs
        v_order(vinf) = 1;
        
        t_in = 1;                   % Order counter
        t0 = 1;                     % Step counter

        fr = 1;                     % Number of candidate ROIs for infection (temp)
        
        % Run dynamics
        while (numel(vinf)>0 && t0<tmax && numel(vinf)<n_con_rois*frate &&...
                numel(fr)>0) 
            % Stop when: 1) no infected nodes; 2) maximum number of simulation
            % steps reached; 3) desired fraction of infected nodes reached; 
            % 4) no more candidate nodes for infection
            
            % Candidate ROIs for infection: 
            f1 = v_neigh_in(vinf,1:kmax)';      % Find neighbours of infected ROIs
            f1 = f1(f1>0); 
            wf1 = v_neigh_w_in(vinf,1:kmax)';   % Probability of each transition
            wf1 = wf1(wf1>0);

            fr = f1(~ismember(f1,vinf));        % Remove already infected ROIs
            wfr = wf1(~ismember(f1,vinf)); 
            nfr = numel(fr);                    % Number of candidates    

            % Select new infected ROIs:
            vc = find(rand(nfr,1)<beta*wfr);    % Idx of new infected nodes
            new_inf = unique(fr(vc)');          % New infected nodes (unique)
            
            % If there are new infected nodes, update stats:
            if (numel(new_inf>0)) 
                t_in = t_in + 1;                % Infection order
                v_order(new_inf) = t_in;        
                v_time(new_inf) = t0;           % Infection time
                vinf = [vinf,new_inf]           % Array of infected ROIs
            end

            t0 = t0+1;

        end
%%  Matrices of infection probability
        v_order(v_order==0) = nr_ROIs+10; %not infected ROIs
        % Infection order:
        order_in = order_in + sparse(v_order,1:nr_ROIs,ones(1,nr_ROIs),...
            nr_ROIs+10,nr_ROIs);
        tt = find(v_time>0);
        xx = 1:nr_ROIs;
        % Infection time:
        time_in = time_in + sparse(v_time(tt),xx(tt),ones(1,numel(tt)),tmax,...
            nr_ROIs);


    end

   
