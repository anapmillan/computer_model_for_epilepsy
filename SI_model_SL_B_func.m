function order_in = SI_model_SL_B_func(bij,seed,nruns)
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
% function order_in = SI_model_SL_B_func(bij,seed,nruns);
% SI_model_B_func runs the SI model in the small beta limit for a binary
% network. It is faster than SI_model_W_func
% INPUTs:
%   bij: binary adjacency matrix
%   seed: seed regions (array)
%   nruns: number of repetitions
%
% OUTPUTs:
%   order_in: probability map of each ROI to become infected at each order 
%             (not normalized).

    nr_ROIs = size(bij,1);
    order_in = zeros(nr_ROIs+10,nr_ROIs);

    ki = sum(bij);
    kmax = max(ki);
    %Create neighbour matrix
    v_neigh = zeros(nr_ROIs,kmax); 
    for i = 1:nr_ROIs
        v_neigh(i,1:ki(i)) = find(bij(i,:));
    end


    % Loop over iterations
    parfor ir = 1:nruns
        v_neigh_in = v_neigh;
        v_order = zeros(nr_ROIs,1); %Matrix of infection times
        vinf = seed;                %Vector of infected ROIs
        v_order(vinf) = 1;
        t = 1;                      %Step counter

        % Candidate ROIs for infection: 
        f1 = v_neigh_in(vinf,:);    %Find neighbours of infected ROIs
        f1 = f1(f1>0);
        fr = f1(~ismember(f1,vinf));%Remove already infected ROIs
        nfr = numel(fr);            %Number of candidates

        % Run dynamics
        while (nfr>0) %Stop when there are no more candidate nodes for infection
            % (all connected nodes are infected)
            t = t+1;
            nup = randi(nfr);       % Select new infected ROI
            vinf = [vinf,fr(nup)];  % Add to infected array
            v_order(fr(nup)) = t;   % Infection time

            % Update vector of candidate ROIs
            f1 = v_neigh(vinf,:);
            f1 = f1(f1>0);
            fr = f1(~ismember(f1,vinf));
            nfr = numel(fr);
        end


    %%  Matrix of infection probability
        v_order(v_order==0) = nr_ROIs+10; %not infected ROIs
        order_in = order_in + ...
            sparse(v_order,1:nr_ROIs,ones(1,nr_ROIs),nr_ROIs+10,nr_ROIs);

    end
