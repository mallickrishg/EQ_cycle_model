function plot_stressperturb_pinnedeqcycles(rcv,V,down,t)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
dipdist=diag(rcv.xc*rcv.dv');
top = 10;
% plot Time and Time-step evolution
subplot(3,2,[1 3])
pcolor(t(1:end-1)./3.15e7,dipdist(1:down:end)/1e3,log10(V(1:end,:)')), shading flat, hold on

box on
set(gca,'YDir','normal');h=colorbar('Location','NorthOutside');
title(h,'log_{10} V'),ylabel('Down-dip distance (km)');
colormap(jet(16))
caxis([-11 -7])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,5);cla;
plot(t(1:end-1)/3.15e7,(max(V(1:end,1:top),[],2)./rcv.Vpl(1)),'m-','LineWidth',2), hold on
xlabel('Time (Yr)','Fontsize',15),ylabel('V/V_{pl}')
axis tight,grid on,ylim([0 2])
set(gca,'FontSize',15,'Color','none')

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(3,2,[2 4]);cla;
imagesc(1:length(t)-1,dipdist(1:down:end)/1e3,log10(V(1:end,:)')), shading flat, hold on

set(gca,'YDir','normal');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('(km)');
caxis([-12 -4])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,6);cla;
plot((1:length(t)-1)/3.15e7,(max(V(1:end,1:top),[],2)./rcv.Vpl(1)),'m-','LineWidth',2), hold on
xlabel('Time Steps','Fontsize',15),ylabel('V/V_{pl}')
axis tight,grid on,ylim([0 2])
set(gca,'FontSize',15,'Color','none')


end