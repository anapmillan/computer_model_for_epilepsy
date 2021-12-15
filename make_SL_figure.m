%% Correlation plot
f_SL = figure;
hh = cell(3,1);

hh{1} = subplot(311);
plot(rho_v, rho_vec_th(:,1),'*--k')
ylabel('PCP')
hh{2} = subplot(312);
plot(rho_v, rho_vec_th(:,3),'*--k')
ylabel('F.Eq.Samp. ROIs')
hh{3} = subplot(313);
plot(rho_v, rho_vec_th(:,5),'*--k')
hold on
plot(rho_v(jmax), rho_vec_th(jmax,5),'^r','markersize',8, 'markerfacecolor','r')
ylabel('Corr. Metric')

for i = 1:3
    xlabel(hh{i}, 'Density')
    if strcmp(scaling_rho,'log')
        set(hh{i},'XScale','log')
    end
end

