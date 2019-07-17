% script to load an appropriate mat file (contains
% rcv,vwparams,perturbparams, downsampling rate, downsampled velocity history, time history)
% Rishav Mallick, EOS, 2019

clear
% inputs
fdir = '~/Documents/results_stressperturb_eqcycles/run';
runnum = 2;
filenum = 17;

% load files
load([fdir num2str(runnum) '/stressperturb_eqcycles_' num2str(filenum) '.mat'])
down = 30;

dipdist = diag(rcv.xc*rcv.dv')./1e3;
%% plot select patch time series as lines
plotdown = 5;

ti = 460; 
tf = 480;

flti = 1;
fltf = 70;

figure(1),clf
% plot_stressperturb_eqcycles(rcv,vwparams,perturbparams,[],Vdown,down,t)
pcolor(t(1:end-1)./3.15e7,dipdist(flti*down:down:fltf*down),log10(Vdown(:,flti:fltf))')
shading flat
colormap(jet(15))
caxis([-11 -8])
box on, grid on
h=colorbar('Location','NorthOutside');
xlabel('Time (yrs)'), ylabel('\zeta_d (km)')
set(gca,'FontSize',15)
xlim([ti tf])

figure(2),clf
vplot = Vdown(:,flti:plotdown:fltf);
nplot = length(vplot(1,:));

for i = 1:nplot
    subplot(nplot,1,i)
    plot(t(1:end-1)./3.15e7,vplot(:,i)./1e-9,'k-','LineWidth',2)
    axis tight, grid on
    ylim([0 2])
    xlim([ti tf])
    title(['\zeta_d = ' num2str(round(dipdist(flti*down + (i-1)*down*plotdown)) ) ' km'])
    set(gca,'FontSize',12)
end

figure(3),clf
colmap = jet(length(vplot(1,:)));
for i = 1:nplot
    plot(t(1:end-1)./3.15e7,vplot(:,i)./1e-9,'-','LineWidth',2,'Color',colmap(i,:)), hold on
    axis tight, grid on
    ylim([0 3])
    xlim([ti tf])
    
end
cb=colorbar;
cb.Label.String = '\zeta_d (km)';
caxis([abs(dipdist(flti*down)) abs(dipdist(down*fltf))])
colormap(colmap)
xlabel('Time (yrs)')
ylabel('V/V_{pl}')
set(gca,'FontSize',15,'YScale','log','Color','none')
ylim([0.01 20])


