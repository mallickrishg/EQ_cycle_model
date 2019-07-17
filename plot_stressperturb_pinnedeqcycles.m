function plot_stressperturb_pinnedeqcycles(rcv,V,down,t)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
dipdist=diag(rcv.xc*rcv.dv');

% plot Time and Time-step evolution
subplot(2,2,[1 3])
pcolor(t(1:end-1)./3.15e7,dipdist(1:down:end)/1e3,log10(V(1:end,:)')), shading flat, hold on

box on
set(gca,'YDir','normal');h=colorbar('Location','NorthOutside');
title(h,'log_{10} V'),ylabel('Down-dip distance (km)');
colormap(jet(16))
caxis([-11 -7])
title(['\Delta \tau = ' num2str(perturbparams.deltau)])
set(gca,'XTickLabel','','FontSize',15)

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(2,2,[2 4]);cla;
imagesc(1:downt:length(t)-1,dipdist(1:down:end)/1e3,log10(V(1:end,:)')), shading flat, hold on

set(gca,'YDir','normal');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('(km)');
caxis([-12 -4])
title(['\Delta \tau = ' num2str(perturbparams.deltau)])
set(gca,'XTickLabel','','FontSize',15)


end