function plot_stressperturb_eqcycles(rcv,vwparams,perturbparams,Y,V,down,t)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
dipdist=diag(rcv.xc*rcv.dv');
tperturb = perturbparams.tperturb;
top = perturbparams.top;
% toppos = floor(top/(rcv.W(1))/(rcv.N/down));
toppos = ceil(top/(rcv.W(1))/(down));

bottom = perturbparams.bot;
botpos = ceil(bottom/(rcv.W(1))/(down));
topvw = vwparams.top;
topvwpos = ceil(topvw/(rcv.W(1))/(down));
botvw = vwparams.bot;
botvwpos = ceil(botvw/(rcv.W(1))/(down));

% down = 1; %spatial downslampling
downt = 1;
% hard-coded Vpl
rcv.Vpl = 1e-9.*ones(rcv.N,1);

% make new colormap
cspec=[cmap('steelblue',100,10,48);...
%     cmap('skyblue',100,45,27);...
    (cmap('skyblue',100,45,27));...
    (cmap('lightgreen',10,59,32));...
    flipud(cmap('orange',100,57,20));...
    flipud(cmap('orangered',100,40,25))];


% plot Time and Time-step evolution
subplot(3,2,[1 3])
pcolor(t(1:end-1)./3.15e7,dipdist(1:down:end)/1e3,log10(V(1:end,:)')), shading flat, hold on
plot([tperturb' tperturb']./3.15e7,get(gca,'YLim'),'k','LineWidth',2)

plot(get(gca,'XLim'),-[top top]./1e3,'m-','Linewidth',2)
plot(get(gca,'XLim'),-[bottom bottom]./1e3,'m-','Linewidth',2)

plot(get(gca,'XLim'),-[topvw topvw]./1e3,'b-','Linewidth',2)
plot(get(gca,'XLim'),-[botvw botvw]./1e3,'b-','Linewidth',2)

box on
set(gca,'YDir','normal');h=colorbar('Location','NorthOutside');
% colormap(cspec);
title(h,'log_{10} V'),ylabel('Down-dip distance (km)');
% caxis([-14 -4])
colormap(jet(16))
caxis([-11 -7])
title(['\Delta \tau = ' num2str(perturbparams.deltau)])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,5);cla;
% plot(t(1:down:end-1)/3.15e7,(mean(V(1:down:end,top:bottom),2)./rcv.Vpl(1)),'m-','LineWidth',2), hold on
plot(t(1:downt:end-1)/3.15e7,(max(V(1:end,toppos:botpos),[],2)./rcv.Vpl(1)),'m-','LineWidth',2), hold on
plot(t(1:downt:end-1)/3.15e7,(max(V(1:end,topvwpos:botvwpos),[],2)./rcv.Vpl(1)),'b-','LineWidth',1)
plot([tperturb' tperturb']./3.15e7,get(gca,'YLim'),'k--','LineWidth',1)
% legend('VS_{perturb}','VW')
% set(legend,'Location','northoutside')
xlabel('Time (Yr)','Fontsize',15),ylabel('V/V_{pl}')
axis tight,grid on,ylim([0 2])
set(gca,'FontSize',15,'Color','none')

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(3,2,[2 4]);cla;
imagesc(1:downt:length(t)-1,dipdist(1:down:end)/1e3,log10(V(1:end,:)')), shading flat, hold on
plot([find(t==tperturb,1) find(t==tperturb,1)],get(gca,'Ylim'),'k','LineWidth',2)

plot(get(gca,'XLim'),-[top top]./1e3,'m-','Linewidth',2)
plot(get(gca,'XLim'),-[bottom bottom]./1e3,'m-','Linewidth',2)

plot(get(gca,'XLim'),-[topvw topvw]./1e3,'b-','Linewidth',2)
plot(get(gca,'XLim'),-[botvw botvw]./1e3,'b-','Linewidth',2)

set(gca,'YDir','normal');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('(km)');
caxis([-12 -4])
title(['\Delta \tau = ' num2str(perturbparams.deltau)])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,6);cla;
%plot(1:downt:length(t)-1,(max(V(1:end,topvwpos:botvwpos),[],2)),'b-','LineWidth',2),hold on
plot(1:downt:length(t)-1,(max(V(1:end,toppos:botpos),[],2)),'m-','LineWidth',2), hold on
xlabel('Time Steps','Fontsize',15),ylabel('V')
axis tight, grid on
plot([find(t==tperturb,1) find(t==tperturb,1)],get(gca,'Ylim'),'k--','LineWidth',1)
set(gca,'FontSize',15,'Color','none','YScale','log')


end