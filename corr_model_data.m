function output = corr_model_data(dd, pattern)
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
% function output = corr_model_data(dd,pattern)
% Measure correlation of spreading pattern (given by dd) with the seizure
% pattern (given by pattern)
% INPUTs:
%   dd:         model data, infection probability map (n_rois x n_orders
%               matrix)
%   pattern:    structure with pattern information
%
% OUTPUTs:
%   output:     cell variable with correlation data: {pcp: pearson's correlation
%               coefficient, frac_act: fraction of sampled rois active in
%               model and data; frac_eq: fraction of sampled rois in the
%               same state (active/inactive) in model and data;
%               n_act_model: fraction of active rois in the model}

    verborrea_tests = 0;
    nrois = size(dd,1);
    nd2 = size(dd,2);

    %% Pattern
    % Sampled ROIs:
    all_samp_rois = unique(pattern.all_erois);
    nrois_sampled = numel(all_samp_rois);

    % Sampled and active ROIs (pattern):
    r_rois = pattern.trois;
    [~,iaux] = unique(r_rois);      % Remove repeated entries
    r_rois = r_rois(sort(iaux));  
    r2_rois = setdiff(all_samp_rois, r_rois); % Sampled and Inactive ROIs in pattern
 
    %% Model
    % Mean activation time
    auxn = sum(dd(:,2:end),2);
    auxn(auxn==0) = 1;
    dd = dd./auxn;
    dm = dd*(0:(nd2-1))';    
    if verborrea_tests>1
            figure
            subplot(121)
            plot(dm)
            subplot(122)
            imagesc(dd)
            pause(1)
    end
    
    %% Overlap pattern - model
    % Vectors of sampled, active, inactive rois (model)
    m_rois = find(dm);                      % Active ROIs (model)
    n_act_model = (numel(m_rois) )/nrois;
    m2_rois = find(~dm);                    % Inactive ROIs in model 
    
    
    mr_rois = intersect(m_rois,r_rois);     % Sampled and Active ROIs (both) 
    mr2_rois = intersect(m2_rois,r2_rois);  % Sampled and Inactive ROIs (both)
    % Fraction of sampled rois that are active in both:
    frac_act = numel(mr_rois)/nrois_sampled;  
    % Fraction of sampled rois that are inactive in both:                                        
    frac_eq = (numel(mr_rois) + numel(mr2_rois)) / nrois_sampled;


    %% Correlation in infection times (only sampled and active ROIs in both)
    % Prepare model data
    dm = dm(mr_rois);   % keep only sampled and active rois
    dmi = unique(dm);   % unique steps for ordering
    i0 = 0;
    a_model = zeros(size(dm)); % array of scaled activation orders (model) 
    for iid = 1:numel(dmi)
        idx = dm==dmi(iid);
        a_model(idx) = i0+mean(1:sum(idx));
        i0 = i0 + sum(idx);
    end

    % Prepare pattern data
    datam = pattern.order;
    ikeep = ismember(r_rois,mr_rois);   % keep only sampled and active rois
    datam = datam(ikeep);
    a_pattern = zeros(size(datam)); % array of scaled activation orders (data)
    datamu = unique(datam);
    i0 = 0;
    for iid = 1:numel(datamu)
        idx = datam==datamu(iid);
        a_pattern(idx) = i0 + mean(1:sum(idx));
        i0 = i0 + sum(idx);
    end


    if numel(a_model) <2 
        pcp = 0; 
    else
        pcp = corr(a_pattern,a_model); 
    end
    
    if verborrea_tests
        figure
        plot(a_pattern)
        hold on
        plot(a_model)
        pause(1)
        title('Re-scaled activation order vectors')
        legend('Pattern','Model')
    end

    output = {pcp, frac_act, frac_eq, n_act_model};
    


end