# computer_model_for_epilepsy

This repository contains a MATLAB library to simulate seizure propagation as a SI spreading process and fit the model to a spreading pattern. The main script is
spreading_master.m

The script allows for two types of SI dynamics and optimization:

1.- Slow Limit (SL): 
  A slow propogation limit that allows for faster calculations. The only free parameter in this case is the network density giving by the network threshold. 
  
2.- Parallel (beta):
  The traditional SI model with a given spreading rate (beta). In this case there are two free parameters, which are fitted in parallel. 
 
 The code includes several options like the possibility of using binary or weighted networks. 
 
 
 The main functions called by spreading_master.m are:
 
    1.- set_data.m (script)
     Specifies the relevant parameters and sets variables, paths, and so on. 
     
    2.- threshold_matrix_rho.m
     Thresholds the full weighted matrix to the desired threshold. 
     
    3.- run_dyn.m
     Takes care of most of the calculations. It calls:
     
      1.- spreading_func.m 
      This selects the adecuate spreading models and runs it:
          1.- SI_model_SL_B_func.m      SL model for a binary network
          2.- SI_model_SL_W_func.m      SL model for a weighted network
          3.- SI_model_beta_W_func.m    Traditional model for a weighted network.
          
      2.- corr_model_data.m
      Calculates the correlation between the seizure pattern and the model, defined as the Pearson's Correlation Coefficient, C.  
     
 The library makes use of three inputs included under "data":
 
 1.- Network structure
 
     Backbone for seizure propagation. An artifitial exemplary network based on the Exponential Distance Rule on a 3D random distribution of centroids (to simulate a brain            network) is included. The size has been set to N=246 in accordance to the Braintome Atlas. The script used to create this network is also included in:
     create_random_EDR_network.m
     
 2.- Epidemic seed
 
     Seed or origin for the epidemics in the model. An artificial exemplary seed is included.
     
 3.- Seizure pattern
 
     Spreading pattern indicating the activation order of the different brain regions during a seizure. This a structure "patttern" including the fields:
     pattern.trois:      (n_sampled_rois x 1) array indicating all sampled ROIs
     pattern.all_erois:  array indicating the indices of sampled active ROIs
     pattern.order:      array indicating the activation order of active sampled
                         ROIs. Size = (n_sampled_active_rois x 1)
     pattern.ntrois:     total number of sampled ROIs
    
    

These programs are distributed by the authors in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use any of these codes, then please cite:
[1]  Ana P Millan et a., "Epidemic models characterize seizure propagation and the effects of epilepsy surgery in individualized brain networks based on MEG and invasive EEG recordings." medRxiv (2021).

[2] Ida Nissen et al. "Optimization of epilepsy surgery through virtual resections on individual structural brain networks." Scientific Reports 11.1 (2021): 1-18.

(c) Ana P. Millan (a.p.millanvidal@amsterdamumc.nl)
