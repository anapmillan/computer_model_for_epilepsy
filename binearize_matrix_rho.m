function bij = binearize_matrix_rho(waij,km0)
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
%   function bij = binearize_matrix_rho(waij,km0)
%   Threshold weighted network, return binary adjacency matrix
%   INPUTs:
%       waij:   full weighted adjacency matrix
%       km0:    target network density
%
%   OUTPUTs:
%       bij:    threshold binary adjancency matrix

    thw = unique(sort(waij));
    thw = [0;thw];
    km = zeros(size(thw));
    for i = 1:numel(thw)
        bij = waij>thw(i);
        bij = (bij+bij') > 0;
        km(i) = mean(sum(bij));            
    end

    [~,ith] = min(abs(km-km0)); 
    bij = waij > thw(ith); 


end