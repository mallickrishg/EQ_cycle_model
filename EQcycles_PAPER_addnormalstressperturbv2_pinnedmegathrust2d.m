%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               EQ CYCLES for a megathrust           %
%               using rate-and-state friction        %
%                                                    %
%               Boundary Integral Method             %
%               Greens'Functions: OKADA 1992         %
% Rishav Mallick, EOS, 2019                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to create long-lived shallow SSE using normal stress perturbations
% (and loading)
% For new SSE paper
% Rishav Mallick, 2019

clear
addpath ~/Dropbox/scripts/unicycle/matlab/
addpath ~/Dropbox/scripts/utils/
import unicycle.*

patchfname='flt_2D_testing_v3.seg';

%% 2. compute greens functions and stress kernels
% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu
G = 30e3; nu = 0.25;
earthModel=greens.okada92(G,nu);
% rcv is the fault
rcv = geometry.receiver(patchfname,earthModel);

% self-stress kernels on the fault
% load stress kernels if they aren't already in the workspace
if ~exist('K','var')
    if exist('stresskernels/Kv3.mat','file')
        disp('Loading K from available file')
        load stresskernels/Kv3.mat
        if length(K)~=rcv.N
            disp('reloading appropriate kernels and fault file')                        
            [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
            save('stresskernels/Kv3.mat','K')
        end
    else
        disp('No K in folder, will create traction kernels')
        [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
        save('stresskernels/Kv3.mat','K')
    end
end
% define dip-distance in (m)
dipdist = -diag(rcv.xc*rcv.dv');
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
fault_width = max(dipdist);
% create velocity weakening region
VW = zeros(rcv.N,1);
topvw = floor(25e3/(fault_width/rcv.N));
botvw = ceil(35e3/(fault_width/rcv.N));
VW(topvw:botvw) = 1;
VW = logical(VW);
% make VW region pinned
rcv.pinnedPosition = VW;
% effective confining pressure on fault (MPa) using depth dependent normal
rcv.sigma = 60.*ones(rcv.N,1);

% frictional parameters 
rcv.a = 1e-2.*ones(rcv.N,1);
% if you want to set the entire thing as VS
rcv.b = rcv.a - 0.5e-2;
rcv.b(VW) = rcv.a(VW) + 0.5e-2;
% static friction coefficient
rcv.mu0 = 0.6*ones(rcv.N,1);
% numer of state parameters
rcv.dgf = 5;
% characteristic weakening distance (m) (use 0.06)
% rcv.l = 0.005.*ones(rcv.N,1);
rcv.l = 0.005.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 1e-9*ones(rcv.N,1);
% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);

% location and rheology of perturbation
% THIS WORKS
% top    = floor(5e3/(fault_width/rcv.N));
% bottom = ceil(25e3/(fault_width/rcv.N));

top    = floor(7e3/(fault_width/rcv.N));
bottom = ceil(18e3/(fault_width/rcv.N));
rcv.b(top:bottom) = rcv.a(top:bottom) - 1e-3;


% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),min([G*rcv.l(top)/rcv.b(top)/rcv.sigma(top) G*rcv.l(topvw-1)/rcv.b(topvw-1)/rcv.sigma(topvw-1)])) 
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(top+1)/rcv.b(top+1))

% add time of perturbation
tperturb = 50.*3.15e7;
trecover_st = 200.*3.15e7;
trecover_end = 225*3.15e7;
% stress perturbation
negstep = 2.1;
deltau = 0;
delsigma = (negstep); % effect of positive delsigma == negative deltau
dsigma_dt = zeros(rcv.N,1);
dsigma_dt(top:bottom) = -negstep./(trecover_end-trecover_st)/1;
% total duration of simulation
tend = 300*3.15e7;

%% Initialize State Vector
Y0=zeros(rcv.N*rcv.dgf,1);
% Fault patches
Y0(1:rcv.dgf:end) = zeros(rcv.N,1);         % slip
Y0(2:rcv.dgf:end) = rcv.sigma.*(rcv.mu0+(rcv.a-rcv.b).*log(rcv.Vpl./rcv.Vo))+G*rcv.Vpl./(2*rcv.Vs); %stress 
Y0(3:rcv.dgf:end) = log(rcv.Vo./rcv.Vpl);      % log(theta Vo / L)
Y0(4:rcv.dgf:end) = log(rcv.Vpl*0.999./rcv.Vo); % log(V/Vo)
Y0(5:rcv.dgf:end) = rcv.sigma;
% initialize the function handle with
% set constitutive parameters
yp=@(t,y) eqcycle_accelerationage_pinnedmegathrust2d_ode(t,y,rcv,K);
% ynp=@(t,y) eqcycle_accelerationage_normload_pinnedmegathrust2d_ode(t,y,rcv,K,dsigma_dt);
ynp=@(t,y) eqcycle_normload_lintaper_2d_ode(t,y,rcv,K,dsigma_dt,trecover_end-trecover_st);

tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-7,'InitialStep',1e-3,'MaxStep',1e7); 
[ti,Yi]=ode45(yp,[0 tperturb],Y0,options);
toc
disp('Stress-Step added')
t = ti;
Y = Yi;
% add a perturbation and continue with solution 
tic
Y0 = Y(end,:)';
th = Y0(3:rcv.dgf:end);
% v+ -> transform to log(v/v0)
sigma2 = rcv.sigma;
sigma2(top:bottom) = sigma2(top:bottom) + delsigma;
sigma1 = rcv.sigma;
Y0((4+(rcv.dgf*(top))):rcv.dgf:(4+(rcv.dgf*(bottom)))) = Y0((4+(rcv.dgf*(top))):rcv.dgf:(4+(rcv.dgf*(bottom)))).*(sigma1(top:bottom)./sigma2(top:bottom)) + ...
    deltau./(rcv.a(top:bottom).*sigma2(top:bottom)) - ...
    (rcv.mu0(top:bottom) + rcv.b(top:bottom).*th(top:bottom)).*(1 - sigma1(top:bottom)./sigma2(top:bottom))./(rcv.a(top:bottom));
Y0(5:rcv.dgf:end) = sigma2;
rcv.sigma(top:bottom) = rcv.sigma(top:bottom) + delsigma;

% evaluate from tperturb to trecover_st
tic
[ti,Yi]=ode45(yp,[0 trecover_st-tperturb],Y0,options);
toc
t = [t;ti+tperturb];
Y = [Y;Yi];

% evaluate from trecover_st to trecover_end
disp('Starting recovery phase')
Y0 = Y(end,:)';
[ti,Yi]=ode45(ynp,[0 trecover_end-trecover_st],Y0,options);
t = [t;ti+trecover_st];
Y = [Y;Yi];
toc
disp('End of recovery phase')

% evaluate from trecover_end to tend
Y0 = Y(end,:)';
[ti,Yi]=ode45(yp,[0 tend-trecover_end],Y0,options);
t = [t;ti+trecover_end];
Y = [Y;Yi];
toc
disp('Finished ODE evaluations')
% Velocities
V=repmat(rcv.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:rcv.dgf:rcv.N*rcv.dgf));
slip = Y(:,1:rcv.dgf:end);
sigmabar = Y(:,5:rcv.dgf:end);
% Maximum Velocity
Vmax=max(V,[],2);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% make new colormap
cpt = readtable('../gmt_figs/eqcycle_continuouscpt.dat');
cspec = [cpt{:,2:4}]./255;
cspec = cspec([1:150],:);

tdown = 1; %spatial downslampling
flti = 10;
fltf = bottom+50;
ti = t(indext_sse)./3.15e7 - 2; 
tf = t(indexend_sse)./3.15e7 - 19;%ti+30;

% PLOT the SSE
Vmax_sse = max(V(find(t==tperturb(1),1):end,top:bottom),[],2);
Vmean_sse = mean(V(find(t==tperturb(1),1):end,top:bottom),2);

indext_sse = find(Vmean_sse > 0.4*rcv.Vpl(1),1) + find(t==tperturb(1),1);
indexend_sse = find((Vmean_sse < 0.4*rcv.Vpl(1)) & (t(find(t==tperturb(1),1):end-1) > t(indext_sse)),1) + find(t==tperturb(1),1);

plotdown = 1;
tdown = 1;
sse_st = t(indext_sse)./3.15e7;
sse_end = t(indexend_sse)./3.15e7;


if isempty(tf)
    tf = max(t)./3.15e7;
end

% flti = 1;
% fltf = bottom+50;%find(abs(dipdist./1e3)>15,1);

figure(10),clf
pcolor(t(1:tdown:(end-1))./3.15e7,dipdist(flti:plotdown:fltf)./1e3,...
    (V(1:tdown:end,flti:plotdown:fltf)./1e-9)')
shading interp
% colormap(jet(30))
colormap(cspec)
caxis([0 2])
hold on
% plot([sse_st sse_st],get(gca,'YLim'),'k--')
% plot([sse_end sse_end],get(gca,'YLim'),'k--')
box on, grid on
xlabel('Time (yrs)'), ylabel('\zeta_d (km)')
set(gca,'FontSize',12,'YDir','reverse')
xlim([ti tf])
% xlim([35 41])
% print('results/Vtransient_focustimeseries','-djpeg','-r600')

figure(11),clf
plot(t./3.15e7,60-mean(sigmabar(:,top:bottom),2),'k-','LineWidth',2)
axis tight
ylim([-2.6 0.2])
yyaxis right
% plot(t./3.15e7,mean(slip(:,top:bottom),2),'k-','LineWidth',2), hold on
plot(t(1:end-1)./3.15e7,mean(V(:,top:bottom),2)./rcv.Vpl(1),'k-','LineWidth',1)
axis tight
% ylim([0.8 2.5])
xlim([ti-3 tf+2])

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
%subplot(3,1,[1 2])
pcolor(t(1:tdown:end-1)./3.15e7,dipdist/1e3,(V(1:tdown:end,:)')./rcv.Vpl(1)), shading flat, hold on

for i = 1:1%length(tperturb)
    plot([tperturb(i) tperturb(i)]./3.15e7,get(gca,'YLim'),'w-','LineWidth',2)
end

plot([ti,tf,tf,ti,ti],dipdist([flti,flti,fltf,fltf,flti])./1e3,'-','Color','w','LineWidth',2)

box on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
% h.Position = [0.1301    0.9010    0.7750    0.0119];
colormap(cspec);
% title(h,'$\frac{V}{V_{pl}}$','interpreter','latex','Fontsize',20),
% colormap(jet(32))
caxis([0 2])
% xlim([tperturb-1*3.15e7 tend]./3.15e7)
xlim([0 tend]./3.15e7)
xlabel('Time (Yr)','Interpreter','latex')
set(gca,'FontSize',15,'YTick',[0,max(dipdist)./4,max(dipdist)/2,3*max(dipdist)./4,max(dipdist)]./1e3,'YTickLabel',{'0';'$\frac{W}{4}$';'$\frac{W}{2}$';'$\frac{3W}{4}$';'W'})
ax = gca;
ax.TickLabelInterpreter = 'latex';


% print('results/Vtransient_timeseries','-djpeg','-r400')
return
%% % Individual patch time series

plotdown2 = 5;

figure(12),clf
vplot = V(:,flti:plotdown2:fltf);
nplot = length(vplot(1,:));
cspec2 = jet(nplot);
for i = 1:nplot
    
    plot(t(1:end-1)./3.15e7,vplot(:,i)./1e-9,'-','LineWidth',2,'Color',cspec2(i,:))
    hold on
    axis tight, grid on
    ylim([0 3])
    xlim([ti tf])
    %title(['\zeta_d = ' num2str(round(dipdist(flti + (i-1)*plotdown2)./1e3) ) ' km'])
end
plot([sse_st sse_st],get(gca,'YLim'),'k--')
plot([sse_end sse_end],get(gca,'YLim'),'k--')

cb=colorbar;
cb.Location = 'northoutside';
cb.Label.String = '\zeta_d (km)';
colormap(cspec2);
caxis([round(dipdist(flti)) round(dipdist(fltf))]./1e3)
xlabel('Time (yrs)')
ylabel('$\frac{V}{V_{pl}}$','Interpreter','latex')
set(gca,'FontSize',12,'Color','none')

%% compute displacements and velocities due to shallow transient ONLY
ndisp = 20;
xplt = linspace(20e3,25e3,ndisp)';
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

figure(14),clf

for i = 1:1:length(dispvelz(:,1))
    
    subplot(2,1,1)
    plot(t(2:end)./3.15e7,(dispz(i,:)'-1*dispz(i,find(t>ti*3.15e7,1)))*1e3,'-','LineWidth',2,'Color',discol(i,:)),hold on
    axis tight, grid on
    xlabel('Time (yr)'),ylabel('Z Displacement (cm)')
    set(gca,'FontSize',12,'Color','none')

%     xlim([ti-0 tf])
    xlim([0 tf])
    
    subplot(2,1,2)
    plot(t(2:end)./3.15e7,dispvelz(i,:)','-','LineWidth',2,'Color',discol(i,:)),hold on
    axis tight, grid on
    xlabel('Time (yr)'),ylabel('Z Velocity (mm/yr)')
    set(gca,'FontSize',12,'Color','none')%'YLim',[-10 10])
    plot([sse_st sse_st],get(gca,'YLim'),'k--')
    plot([sse_end sse_end],get(gca,'YLim'),'k--')
%     xlim([ti-0 tf])
    xlim([0 tf])
    ylim([-15 0])
    


end
return
%% plot select velocity profiles
tplot1 = 120*3.15e7;
tplot2 = 210*3.15e7;

figure(15),clf
% set(gcf,'Renderer','painters')
tindex = find(abs(t-tplot1) == min(abs(t-tplot1)));

% % figure(16),clf
% % set(gcf,'Renderer','painters')
% % plot(dipdist,V(tindex+1,:)./rcv.Vpl(1),'-','Color',rgb('gray'),'LineWidth',2)
% axis tight
% ylim([0 1.5])
% % set(gca,'FontSize',12,'XTick',[0,max(dipdist)./2,max(dipdist)],'XTickLabel',{'0';'$\frac{W}{2}$';'W'},'YTick',[0:0.5:1.5])
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% 
% ylabel('$\left(\frac{V}{V_{pl}}\right)$','Interpreter','latex','FontSize',20)
% print('results/Vtransientprofile','-depsc','-r300')
plot(V(tindex+1,:)./rcv.Vpl(1),dipdist,'-','Color',rgb('gray'),'LineWidth',2), hold on
axis tight
xlim([0 1.5])
set(gca,'FontSize',12,'YTick',[0:0.25:1].*max(dipdist),'YTickLabel','',...
    'XTick',[0:0.5:1.5],'XTickLabel','','YDir','reverse')
ax = gca;
ax.TickLabelInterpreter = 'latex';
% xlabel('$\left(\frac{V}{V_{pl}}\right)$','Interpreter','latex','FontSize',20)

tindex = find(abs(t-tplot2) == min(abs(t-tplot2)));
plot(V(tindex+1,:)./rcv.Vpl(1),dipdist,'-','Color',rgb('red'),'LineWidth',2), hold on

% make figure for a-b
ax1=gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(rcv.a-rcv.b,dipdist,'Parent',ax2,'LineWidth',2,'Color',rgb('royalblue'))
hold on
% line(-5*(rcv.a-rcv.b),dipdist,'Parent',ax2,'LineWidth',2,'Color',rgb('royalblue'))
line(rcv.sigma.*7e-4,dipdist,'Linewidth',2,'Color',rgb('forestgreen'))
set(gca,'YDir','reverse','YTickLabel','','XTickLabel','','XTick','')
% print('results/Vlockedprofile','-depsc','-r300')
% 