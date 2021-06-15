%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               EQ CYCLES for a megathrust           %
%               using rate-and-state friction        %
%                                                    %
%               Boundary Integral Method             %
%               Greens'Functions: OKADA 1992         %
% Rishav Mallick, EOS, 2019                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kinematically impose normal stress increase (fast expulsion of fluids)
% and return (slow migration of fluids)
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
topvw = floor(10e3/(fault_width/rcv.N));
botvw = ceil(25e3/(fault_width/rcv.N));
% VW(topvw:botvw) = 1;
VW = logical(VW);
% make VW region pinned
rcv.pinnedPosition = VW;
% effective confining pressure on fault (MPa) using depth dependent normal
rcv.sigma = 30.*ones(rcv.N,1);

% frictional parameters 
rcv.a = 1e-2.*ones(rcv.N,1);
% if you want to set the entire thing as VS
rcv.b = rcv.a - 5e-3;

% static friction coefficient
rcv.mu0 = 0.6*ones(rcv.N,1);
% numer of state parameters
rcv.dgf = 5;
% characteristic weakening distance (m) (use 0.06)
rcv.l = 0.02.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 1e-9*ones(rcv.N,1);
% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);

% location and rheology of perturbation
top    = floor(30e3/(fault_width/rcv.N));
bottom = ceil(40e3/(fault_width/rcv.N));
rcv.b(top:bottom) = rcv.a(top:bottom) - 3e-3;

% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),min([G*rcv.l(top)/rcv.b(top)/rcv.sigma(top) G*rcv.l(topvw-1)/rcv.b(topvw-1)/rcv.sigma(topvw-1)])) 
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(top+1)/rcv.b(top+1))

% stress perturbation
tperturb = 20; %start first perturbation
Trec = 10; %recurrence time of SSE in years
Ncycles = 7;

deltau = 0; % traction change in MPa
delsigma = 2; % (MPa) effect of positive delsigma == negative deltau

% normal stress loading (fluid return over Trec)
dsigma_dt = zeros(rcv.N,1);
dsigma_dt(top:bottom) = -delsigma./(Trec*3.15e7);


% total duration of simulation
tend = 101*3.15e7;


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
ynp=@(t,y) eqcycle_accelerationage_normload_pinnedmegathrust2d_ode(t,y,rcv,K,dsigma_dt);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-7,'InitialStep',1e-3,'MaxStep',1e7); 
[ti,Yi]=ode45(yp,[0 tperturb*3.15e7],Y0,options);
toc
disp('Stress-Step added')
t = ti;
Y = Yi;
%% add a perturbation and continue with solution (loop)
tic
for ii = 1:Ncycles
    Y0 = Y(end,:)';
    th = Y0(3:rcv.dgf:end);
    % compute new normal stresses
    sigma2 = Y0(5:rcv.dgf:end);
    sigma2(top:bottom) = sigma2(top:bottom) + delsigma;
    sigma1 = Y0(5:rcv.dgf:end);
    rcv.sigma = sigma2;
    % v+ from v- and change in sigma
    %Y0(4*top:rcv.dgf:bottom*rcv.dgf)
    vminus = Y0(4:rcv.dgf:end);
    vplus = vminus;
    vplus(top:bottom) = vminus(top:bottom).*(sigma1(top:bottom)./sigma2(top:bottom)) + ...
        deltau./(rcv.a(top:bottom).*sigma2(top:bottom)) - ...
        (rcv.mu0(top:bottom) + rcv.b(top:bottom).*th(top:bottom)).*(1 - sigma1(top:bottom)./sigma2(top:bottom))./(rcv.a(top:bottom));
    Y0(4:rcv.dgf:end) = vplus;
    
    Y0(5:rcv.dgf:end) = sigma2;
    
    % evaluate from tperturb to trecover_st
    tic
    [ti,Yi]=ode45(ynp,[0 Trec*3.15e7],Y0,options);
    toc
    t = [t;ti+t(end)];
    Y = [Y;Yi];
end

Y0 = Y(end,:)';
[ti,Yi]=ode45(yp,[0 tend-t(end)],Y0,options);
t = [t;ti+t(end)];
Y = [Y;Yi];
toc
disp('End of Simulation')

% Extract important quantities
V=repmat(rcv.Vo',size(Y,1),1).*exp(Y(1:end,4:rcv.dgf:rcv.N*rcv.dgf));
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

tdown = 1; %temporal downslampling

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
%subplot(3,1,[1 2])
pcolor(t(1:tdown:end)./3.15e7,dipdist/1e3,(V(1:tdown:end,:)')./rcv.Vpl(1)), shading flat, hold on

% for i = 1:1%length(tperturb)
%     plot([tperturb(i) tperturb(i)]./3.15e7,get(gca,'YLim'),'w--','LineWidth',1)
% end
plot(get(gca,'XLim'),[dipdist(top) dipdist(top)]./1e3,'w-','Linewidth',.5)
plot(get(gca,'XLim'),[dipdist(bottom) dipdist(bottom)]./1e3,'w-','Linewidth',.5)
box on
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
h.TickLabelInterpreter = 'latex';
colormap(cspec);
title(h,'$\frac{V}{V_{pl}}$','interpreter','latex','Fontsize',20),
% colormap(jet(32))
caxis([0 2])
xlabel('Time (Yr)','Interpreter','latex')
set(gca,'FontSize',15,'YTick',[0,max(dipdist)./4,max(dipdist)/2,3*max(dipdist)./4,max(dipdist)]./1e3,'YTickLabel',{'0';'$\frac{W}{4}$';'$\frac{W}{2}$';'$\frac{3W}{4}$';'W'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
ylabel('$\frac{\zeta_d}{W}$','Interpreter','latex','FontSize',20);

figure(2),clf
plot(t./3.15e7,max(V./rcv.Vpl(1),[],2),'LineWidth',2)
xlabel('Time (Yr)','Interpreter','latex')
ylabel('$\frac{V}{V_{pl}}$','Interpreter','latex')
set(gca,'FontSize',15)
ax = gca;
ax.TickLabelInterpreter = 'latex';
axis tight
grid on

yyaxis right
plot(t./3.15e7,max(sigmabar(:,top:bottom),[],2),'LineWidth',1)
ylabel('$\bar{\sigma}$ (MPa)','Interpreter','latex')
return
%% normal stress and slip history
figure(3)
subplot(2,1,1)
pcolor(t(1:tdown:end)./3.15e7,dipdist/1e3,(sigmabar(1:tdown:end,:)')), shading flat, hold on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
caxis([max(rcv.sigma)-delsigma max(rcv.sigma)+delsigma])

subplot(2,1,2)
pcolor(t(1:tdown:end)./3.15e7,dipdist/1e3,(slip(1:tdown:end,:)')), shading flat, hold on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
colormap jet