% script to load an appropriate mat file (contains
% rcv,vwparams,perturbparams, downsampling rate, downsampled velocity history, time history)
% must contain the undownweighted Y,V,rcv,t etc
% Rishav Mallick, EOS, 2019
%% PLOT the SSE
Vmax_sse = max(V(find(t==tperturb(1),1):end,top:bottom),[],2);
Vmean_sse = mean(V(find(t==tperturb(1),1):end,top:bottom),2);

% indext_sse = find(Vmax_sse > 0.3*rcv.Vpl(1),1) + find(t==tperturb(1),1);
indext_sse = find(Vmean_sse > 0.3*rcv.Vpl(1),1) + find(t==tperturb(1),1);

% indexend_sse = find((Vmax_sse < 0.3*rcv.Vpl(1)) & (t(find(t==tperturb(1),1):end-1)>t(indext_sse)),1) + find(t==tperturb(1),1);
indexend_sse = find((Vmean_sse < 0.3*rcv.Vpl(1)) & (t(find(t==tperturb(1),1):end-1)>t(indext_sse)),1) + find(t==tperturb(1),1);

plotdown = 5;
tdown = 3;
% vplot = Y(2:end,2:rcv.dgf:end);

ti = t(indext_sse)./3.15e7 - 5; 
tf = t(indexend_sse)./3.15e7;%ti+30;

flti = 1;
fltf = find(abs(dipdist./1e3)>30,1);

figure(10),clf

% plot_stressperturb_eqcycles(rcv,vwparams,perturbparams,[],Vdown,down,t)
% pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,log10(V(1:tdown:end,flti:plotdown:fltf))')
pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,(V(1:tdown:end,flti:plotdown:fltf)./1e-9)')
% pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,(vplot(1:tdown:end,flti:plotdown:fltf)') - vplot(find(t>ti*3.15e7,1),flti:plotdown:fltf)')

% caxis([-1 1]*.3)
shading interp
colormap(jet(30))
% colormap bluewhitered
% caxis([-11 -7])
caxis([0 2])


box on, grid on
h=colorbar('Location','NorthOutside');
xlabel('Time (yrs)'), ylabel('\zeta_d (km)')
set(gca,'FontSize',15,'YDir','reverse')
xlim([ti tf])
% xlim([35 41])
%% compute displacements and velocities due to shallow transient ONLY
ndisp = 20;
xplt = linspace(7e3,20e3,ndisp)';
[~,Gd]=rcv.displacementKernels([xplt,0.*xplt,0.*xplt],3);
Gdz = Gd(3:3:end,:);
discol = jet(ndisp);

slip_t0 = (Y(2:end,1:rcv.dgf:end)');
backslip = 0*rcv.Vpl*(t(2:end)');
% slip_t0(find(abs(dipdist./1e3)>30,1):end,:) = 0; % remove displacements/velocities from VW region

V_t0 = V';
% V_t0(find(abs(dipdist./1e3)>30,1):end,:) = 0; % remove displacements/velocities from VW region

dispz = Gdz*slip_t0 - 1*Gdz*backslip;
dispvelz = 3.15e10.*(Gdz*(V_t0-rcv.Vpl));
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

scf = 3e1;


figure(11),clf
surf(Xp,Tp',Zp,cp), hold on
alpha 0.5
axis tight, shading flat, box on
caxis([0 2])
colormap(jet(30))

ylim([ti tf])

% xlim(minmax(xp))

% view(40,25)
view(90,90)

% daspect([1 1 1])
xlabel('X (km)'),zlabel('Z (km)'),ylabel('T (yrs)')
h=colorbar('Location','EastOutside');
h.Label.String = 'V/V_{pl}';
set(gca,'FontSize',15,'Color','none')

% plot scaled GPS time series as wiggles on map
plot3(xplt./1e3,ti.*ones(size(xplt)),zeros(size(xplt)),'kv','MarkerFaceColor','k')
for i = 1:1:length(xplt)%10:2:30
    dispts = dispz(i,:)'-1*dispz(i,find(t>ti*3.15e7,1));
    
    plot3(xplt(i)/1e3.*ones(length(t(2:end)),1) - dispts.*scf,t(2:end)./3.15e7,zeros(length(t(2:end)),1),'-','Color',discol(i,:),'LineWidth',2)
end
xlim([min(xp) xplt(i)./1e3 + 3])

% print('Figures/3d_transient_timeseries','-djpeg','-r300')

%% Individual patch time series

plotdown2 = 15;

figure(12),clf
vplot = V(:,flti:plotdown2:fltf/2);
nplot = length(vplot(1,:));

for i = 1:nplot
    subplot(nplot,1,i)
    plot(t(1:end-1)./3.15e7,vplot(:,i)./1e-9,'k-','LineWidth',2)
    axis tight, grid on
    ylim([0 10])
    xlim([ti tf])
    title(['\zeta_d = ' num2str(round(dipdist(flti + (i-1)*plotdown2)./1e3) ) ' km'])
    set(gca,'FontSize',12)
end



%% displacement time series (on surface)
figure(14),clf
% ti = 9;
for i = 1:4:length(dispvelz(:,1))
    
    subplot(2,1,1)
    plot(t(2:end)./3.15e7,(dispz(i,:)'-1*dispz(i,find(t>ti*3.15e7,1)))*1e3,'-','LineWidth',2,'Color',discol(i,:)),hold on
    axis tight, grid on
    xlabel('Time (yr)'),ylabel('Z Displacement (cm)')
    set(gca,'FontSize',12,'Color','none')

    xlim([ti-10 tf])
    
    subplot(2,1,2)
    plot(t(2:end)./3.15e7,dispvelz(i,:)','-','LineWidth',2,'Color',discol(i,:)),hold on
    axis tight, grid on
    xlabel('Time (yr)'),ylabel('Z Velocity (mm/yr)')
    set(gca,'FontSize',12,'Color','none')%'YLim',[-10 10])
    plot([sse_st sse_st],get(gca,'YLim'),'k--')
    plot([sse_end sse_end],get(gca,'YLim'),'k--')
    xlim([ti-100 tf])
    ylim([-1 0].*15)
    


end

%% 3-d view of GPS vertical timeseries
