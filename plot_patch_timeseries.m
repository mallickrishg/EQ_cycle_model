% script to load an appropriate mat file (contains
% rcv,vwparams,perturbparams, downsampling rate, downsampled velocity history, time history)
% Rishav Mallick, EOS, 2019

clear
% inputs
fdir = '~/Documents/results_stressperturb_eqcycles/run';
runnum = 23;
filenum = 2;

% load files
load([fdir num2str(runnum) '/stressperturb_eqcycles_' num2str(filenum) '.mat'])
down = 30;

dipdist = diag(rcv.xc*rcv.dv')./1e3;
%% plot select patch time series as lines
plotdown = 4;
colvals = [];
ti = 900; 
tf = 960;

flti = 1;
fltf = 40;

figure(1),clf
% plot_stressperturb_eqcycles(rcv,vwparams,perturbparams,[],Vdown,down,t)
subplot(2,1,1)
pcolor(t(1:end-1)./3.15e7,dipdist(flti*down:down:fltf*down),log10(Vdown(:,flti:fltf))')
shading flat
colormap(jet(16))
caxis([-10 -8])
box on, grid on
h=colorbar('Location','NorthOutside');
h.Label.String = 'log_{10} V';
xlabel('Time (yrs)'), ylabel('\zeta_d (km)')
set(gca,'FontSize',15)
xlim([ti tf])

subplot(2,1,2)
pcolor(t(1:end-1)./3.15e7,dipdist(flti*down:down:fltf*down),(Vdown(:,flti:fltf))'./rcv.Vpl(1)), shading flat
colormap(jet(16))
caxis([0 2])
box on, grid on
h=colorbar('Location','NorthOutside');
h.Label.String = 'V/V_{pl}';
xlabel('Time (yrs)'), ylabel('\zeta_d (km)')
set(gca,'FontSize',15)
xlim([ti tf])



% time series
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
    colvals(i) = dipdist(flti*down + (i-1)*down*plotdown);
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
caxis(minmax(abs(colvals)) + [-1 1])
colormap(colmap)
cb.Ticks = round(abs(colvals));
xlabel('Time (yrs)')
ylabel('V/V_{pl}')
set(gca,'FontSize',15,'YScale','log','Color','none')
ylim([0.01 20])


