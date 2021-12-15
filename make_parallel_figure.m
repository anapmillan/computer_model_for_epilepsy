f_parallel = figure;

hh{1} = subplot(131);
imagesc(squeeze(rho_vec_parallel(:,:,1)))
title('PCP')

hh{2} = subplot(132);
imagesc(squeeze(rho_vec_parallel(:,:,3)))
title('Frac. Eq. ROIs')

hh{3} = subplot(133);
imagesc(squeeze(rho_vec_parallel(:,:,5)))
title('Correlation metric')

for i = 1:3
    ylabel(hh{i},'Threshold')
    xlabel(hh{i},'Spreading rate')
    colorbar(hh{i})
    set(hh{i},'YTick',1:numel(rho_v), 'YTickLabel',num2str(rho_v(1:end)','%.4f'));
    set(hh{i},'XTick',1:numel(rho_v), 'XTickLabel',num2str(beta_v(1:end)','%.4f'));
    set(hh{i},'YTickLabelRotation',80)
end