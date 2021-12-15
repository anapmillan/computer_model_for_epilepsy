function order_in = SI_model_SL_W_func(aij,seed,nruns)
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
% function order_in = SI_model_SL_W_func(aij, seed, nruns);
% SI_model_W_func runs the SI model in the small beta limit for a weighted
%   network. 
% INPUTs:
%   aij: weighted adjacency matrix (nrois x nrois matrix)
%   seed: seed regions (1xNR array)
%   nruns: number of repetitions
%
% OUTPUTs:
%   order_in: probability map of each ROI to become infected at each step (not normalized).

    nr_ROIs = size(aij,1);              % Number of ROIs 

    order_in = zeros(nr_ROIs+10,nr_ROIs);

    % Create Degree array and neighbour matrix
    bij = aij>0;                        % Binary matrix
    ki = sum(bij);                      % Degree array
    kmax = max(ki);
    v_neigh = zeros(nr_ROIs,kmax);      % Neighbour matrix (from binary network)
    v_neigh_w = zeros(nr_ROIs,kmax);    % Weighted neihbour matrix (for probabilities)
    for i = 1:nr_ROIs
        v_neigh(i,1:ki(i)) = find(bij(i,:));
        v_neigh_w(i,1:ki(i)) = aij(i,v_neigh(i,1:ki(i)));
    end

    % Loop over iterations
    parfor ir = 1:nruns    
        v_neigh_in = v_neigh;
        v_neigh_w_in = v_neigh_w;
        
        v_order = zeros(nr_ROIs,1); %Matrix of infection times
        vinf = seed';                %Vector of infected ROIs
        v_order(vinf) = 1;
        t = 1;                      %Step counter
        
        % Candidate ROIs for infection: 
        f1 = v_neigh_in(vinf,:);    %Find neighbours of infected ROIs
        f1 = f1(f1>0);
        wf1 = v_neigh_w_in(vinf,:); %Probability of each transition
        wf1 = wf1(wf1>0);

        fr = f1(~ismember(f1,vinf)); %Remove already infected ROIs
        wfr = wf1(~ismember(f1,vinf));
        nfr = numel(fr);

        % Run dynamics
        while (nfr>0) %Stop when there are no more candidate nodes for infection
            % (all connected nodes are infected)
            t = t+1;   
            % Select new infected ROI
            auxn = sum(wfr) * rand; 
            % sum(wfr): total transition rate
            % Find ROI corresponding to auxn:
            dx = 0;
            ix = 1;
            while dx < auxn
                dx = dx + wfr(ix);
                ix = ix + 1;
            end
            nup = ix - 1;           % Index of newly infected ROI

            % New infected ROI:
            vinf = [vinf,fr(nup)];
            v_order(fr(nup)) = t;

            % Update vector of candidate ROIs            
            f1 = v_neigh(vinf,:);
            f1 = f1(f1>0);
            fr = f1(~ismember(f1,vinf));
            wf1 = v_neigh_w(vinf,:);
            wf1 = wf1(wf1>0);
            wfr = wf1(~ismember(f1,vinf));
            nfr = numel(fr);
        end


    %%  Matrix of infection probability
        v_order(v_order==0) = nr_ROIs+10; %not infected ROIs
        order_in = order_in + sparse(v_order,...
            1:nr_ROIs,ones(1,nr_ROIs),nr_ROIs+10,nr_ROIs);
    end