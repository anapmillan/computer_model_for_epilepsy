function [tot_pob, m_mass, tot_pob_time, m_mass_t] = ...
    spreading_func(aij, seed, nrunss, str_name, pout, ...
    w_conn, w_dyn, beta, tmax, save_data, mfig)
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
%function [tot_pob, m_mass, pov3, m_mass_t] = ...
%    spreading_func(aij, seed, nrunss, str_name, pout, ...
%    w_conn, w_dyn, beta, tmax, mfig)
% INPUTS:
%   aij: adj or weight matrix
%   seed: seed ROIs
%   nrunss: number of runs of epidemic model (for averaging)
%   str_name, pout: name tag and directory for filemanes
%   w_conn, w_dyn: conn type (B/W) and model type (SI, beta)
%   beta, tmax, gamma: parameters for model
%   save_data : switch to turn on/off saving data in file (default: 1)
%   mfig: switch to get some plots of spreading (detault: 0)
%
% OUTPUTS:
%   tot_pob: probability map of a region getting infected at a given order
%   m_mass: fraction of infected rois at each infection order
%   tot_pob_time: probability map of a region being infected at a given
%      step (only beta)
%   m_mass_t: fraction of infected rois at each step (only beta)


    do_verborrea = 0;

    time_out = 0;
    order_out = 0;
    tot_pob_time = 0;
    m_mass_t = 0;
    m_mass = 0;
    
    if nargin<10
        save_data = 1;
        if nargin<11
            mfig = 0;
        end
    end
    
    %% Recover Some parameters
    nrois = size(aij,1);
    rkm = sum(aij(:)>0)/nrois;
 
    str_name_path = sprintf('%s/dynamics_%s',pout,str_name);
  
    
    %% Select model and RUN
    if strcmp(w_dyn,'SI') % Model 1  
        if do_verborrea; fprintf('Beta -> 0 limit\n'); end      
        if strcmp(w_conn,'W') %Weighted matrix
            order_in = SI_model_SL_W_func(aij,seed,nrunss);
        elseif strcmp(w_conn,'B') %Binary matrix
            order_in = SI_model_SL_B_func(aij,seed,nrunss);
        end
        tot_pob = order_in(1:nrois,:)';
        
    elseif strcmp(w_dyn,'beta') %Model 2
        if do_verborrea; fprintf('Finite beta mode\n'); end
        % Two outputs: activation order and time
        [order_in,time_in] = SI_model_beta_W_func(...
            aij,seed,beta,nrunss,tmax);
        tot_pob = order_in';
        tot_pob_time = time_in';       
    end
    
    tot_pob = tot_pob/nrunss;
    tot_pob_time = tot_pob_time/nrunss;
    
    
    
    %% Output files & Write
    if save_data
        if strcmp(w_dyn,'SI')
            n_out = sprintf('%s_km%.3f_nruns%d',str_name_path,rkm, nrunss);        
        elseif strcmp(w_dyn,'beta')
            n_out = sprintf('%s_km%.3f_beta%.4f_tmax%d_nruns%d',...
                str_name_path,rkm, beta, tmax, nrunss);
        elseif strcmp(w_dyn,'SIR')
            n_out = sprintf('%s_km%.3f_beta%.4f_gamma%.4f_tmax%d_nruns%d',...
                str_name_path,rkm, beta, gamma, tmax, nrunss);
        else
            err('Incorrect model %s\n', w_dyn)
        end
        fprintf('Saving spreading data into %s\n', n_out)
        dlmwrite(sprintf('%s.txt',n_out),tot_pob,' ');  

    
        %% Some plots

        if mfig
            % 2D Probability map
            aux_color = 0.2;
            tot_pob_time = tot_pob;
            tot_pob_time(seed,1)=0;
            xmax = max(find(sum(tot_pob)>0));
            %xmax = 10 * ceil(xmax/10);
            pmax = max(max(tot_pob_time));

            figure
            imagesc(tot_pob)
            caxis([0 0.2*pmax])
            colorbar
            axis([0 xmax 0 nrois]);
            xlabel('Activation time');
            ylabel('ROIs');
            title(sprintf('%s %s',str_name,rkm));
            print('-dpng',sprintf('%s_unsorted.png',n_out));

            %% Spreading curve
            m_mass_in = mean(order_in,2)/nrunss;
            m_mass = zeros(size(m_mass_in));
            for i = 1:(numel(m_mass))
                m_mass(i) = sum(m_mass_in(1:i));
            end

            figure
            if strcmp(w_dyn,'beta')
                m_mass_t_in = mean(time_in,2)/nrunss;
                tmax_r = max(find(m_mass_t_in>0))+1;
                m_mass_t = zeros(tmax_r,1);
                for i = 1:tmax_r
                    m_mass_t(i) = sum(m_mass_t_in(1:i));
                end
                subplot(211)
                plot(m_mass_t, 'linewidth', 2)
                subplot(212)
            end
            plot(m_mass, 'linewidth', 2)
            n_mass_plot= sprintf('%s_mass.png', n_out);
            print('-dpng',n_mass_plot);

        end
    end

end