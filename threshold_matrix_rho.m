function aij = threshold_matrix_rho(wij,rho,N)
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
%   function aij = threshold_matrix_rho(wij,rho,N,is_eu)
%   Threshold weight matrix to desired threshold
%   INPUTS:
%       wij:    full weight matrix
%       rho:    threshold value
%       N:      network size (number of nodes)
%
%   OUTPUTs:
%       aij:    thresholded weighted matrix
    kmo = rho*N;
    bij = binearize_matrix_rho(wij,kmo);
    bij = (bij + bij')>0;
    aij = wij .* bij;
