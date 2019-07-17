% script to load an appropriate mat file (contains
% rcv,vwparams,perturbparams, downsampling rate, downsampled velocity history, time history)
% must contain the undownweighted Y,V,rcv,t etc
% Rishav Mallick, EOS, 2019

%% plot select patch time series as lines
plotdown = 3;
tdown = 5;
vplot = Y(2:end,2:rcv.dgf:end);

ti = 467; 
tf = 475;

flti = 1;
fltf = 1000;
%%
figure(1),clf

% plot_stressperturb_eqcycles(rcv,vwparams,perturbparams,[],Vdown,down,t)
% pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,log10(V(1:tdown:end,flti:plotdown:fltf))')
% pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,(V(1:tdown:end,flti:plotdown:fltf)./1e-9)')
pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,(vplot(1:tdown:end,flti:plotdown:fltf)') - vplot(find(t>ti*3.15e7,1),flti:plotdown:fltf)')

caxis([-1 1]*.3)
shading interp
% colormap(jet(15))
colormap bluewhitered
% caxis([-11 -7])
% caxis([0 5])


box on, grid on
h=colorbar('Location','NorthOutside');
xlabel('Time (yrs)'), ylabel('\zeta_d (km)')
set(gca,'FontSize',15)
xlim([ti tf])

%% plot 3-D view
% plotdown = 10;
% tdown = 20;

tp = t(1:tdown:(end-1))./3.15e7;
xp = rcv.xc(flti:plotdown:fltf,1)./1e3;
zp = rcv.xc(flti:plotdown:fltf,3)./1e3;
cp = (V(1:tdown:end,flti:plotdown:fltf)./1e-9)';

nt = length(tp);
nf = length(xp);

Tp = repmat(tp,1,nf);
Xp = repmat(xp,1,nt);
Zp = repmat(zp,1,nt);

figure(11),clf
surf(Xp,Tp',Zp,cp)

axis tight, shading interp
caxis([0 3])
colormap(jet(30))
ylim([ti tf])
view(40,25)
daspect([1 .3 1])
xlabel('X (km)'),zlabel('Z (km)'),ylabel('T (yrs)')
h=colorbar('Location','EastOutside');
h.Label.String = 'V/V_{pl}';
set(gca,'FontSize',15,'Color','none')

print('Figures/3d_transient_timeseries','-djpeg','-r300')

%% Individual patch time series

plotdown2 = 100;

figure(12),clf
vplot = V(:,flti:plotdown2:fltf);
nplot = length(vplot(1,:));

for i = 1:nplot
    subplot(nplot,1,i)
    plot(t(1:end-1)./3.15e7,vplot(:,i)./1e-9,'k-','LineWidth',2)
    axis tight, grid on
    ylim([0 3])
    xlim([ti tf])
    title(['\zeta_d = ' num2str(round(dipdist(flti + (i-1)*plotdown2)./1e3) ) ' km'])
    set(gca,'FontSize',12)
end

%% compute displacements and velocities due to shallow transient ONLY
ndisp = 50;
xplt = linspace(0.1e3,120e3,ndisp)';
[~,Gd]=rcv.displacementKernels([xplt,0.*xplt,0.*xplt],3);
Gdz = Gd(3:3:end,:);
discol = hot(ndisp);

slip_t0 = (Y(2:end,1:rcv.dgf:end)');
slip_t0(1000:end,:) = 0; % remove displacements/velocities from VW region
V_t0 = V';
V_t0(1000:end,:) = 0; % remove displacements/velocities from VW region

dispz = Gdz*slip_t0;
dispvelz = 3.15e10.*(Gdz*V_t0);
%%
figure(4),clf
for i = 13:1:30%length(dispvelz(:,1))
    
    subplot(1,1,1)
    plot(t(2:end)./3.15e7,dispz(i,:)'-dispz(i,find(t>ti*3.15e7,1)),'-','LineWidth',2,'Color',discol(i,:)),hold on
    axis tight, grid on
    xlabel('Time (yr)'),ylabel('Z Displacement (m)')
    set(gca,'FontSize',15,'Color','none')
    xlim([ti tf])
    
%     subplot(2,1,2)
%     plot(t(2:end)./3.15e7,dispvelz(i,:)','-','LineWidth',2,'Color',discol(i,:)),hold on
%     axis tight, grid on
%     xlabel('Time (yr)'),ylabel('Z Velocity (mm/yr)')
%     set(gca,'FontSize',15,'Color','none','YLim',[-10 10])
%     xlim([ti tf])

end

%% 3-d view of GPS vertical timeseries
