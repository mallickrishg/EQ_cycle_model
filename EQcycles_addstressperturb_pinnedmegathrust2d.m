%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               EQ CYCLES for a megathrust           %
%               using rate-and-state friction        %
%                                                    %
%               Boundary Integral Method             %
%               Greens'Functions: OKADA 1992         %
% Rishav Mallick, EOS, 2019                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

patchfname='flt_2D_testing.seg';

patchdeepname='flt_2D_deep.seg';

%% 2. compute greens functions and stress kernels
% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu
G = 30e3; nu = 0.25;
earthModel=greens.okada92(G,nu);
% rcv is the fault
rcv = geometry.receiver(patchfname,earthModel);
rcvdeep = geometry.receiver(patchdeepname,earthModel);
% self-stress kernels on the fault
% load stress kernels if they aren't already in the workspace
if ~exist('K','var')
    if exist('stresskernels/K.mat','file')
        disp('Loading K from available file')
        load stresskernels/K.mat
        if length(K)~=rcv.N
            disp('reloading appropriate kernels and fault file')                        
            [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
            save('stresskernels/K.mat','K')
        end
    else
        disp('No K in folder, will create traction kernels')
        [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
        save('stresskernels/K.mat','K')
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
topvw = floor(50e3/(fault_width/rcv.N));
botvw = ceil(80e3/(fault_width/rcv.N));
VW(topvw:botvw) = 1;
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
rcv.mu0 = 0.1*ones(rcv.N,1);
% numer of state parameters
rcv.dgf=4;
% characteristic weakening distance (m)
rcv.l = 0.03.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 1e-9*ones(rcv.N,1);
% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);

% location and rheology of perturbation
top    = floor(1e3/(fault_width/rcv.N));
bottom = ceil(20e3/(fault_width/rcv.N));
rcv.b(top:bottom) = rcv.a(top:bottom) - 0.1e-3;


% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),min([G*rcv.l(top)/rcv.b(top)/rcv.sigma(top) G*rcv.l(topvw-1)/rcv.b(topvw-1)/rcv.sigma(topvw-1)])) 
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(1)/rcv.b(1))
fprintf(1,'a/b in VW regions = %.2f\n',mean(rcv.a(VW)./rcv.b(VW)))

% add time of perturbation
tperturb = 100*3.15e7;
% stress perturbation
deltau = [-0.5];
% total duration of simulation
tend = 300*3.15e7;
%% Initialize State Vector
Y0=zeros(rcv.N*rcv.dgf,1);
% Fault patches
Y0(1:rcv.dgf:end) = zeros(rcv.N,1);         % slip
Y0(2:rcv.dgf:end) = rcv.sigma.*(rcv.mu0+(rcv.a-rcv.b).*log(rcv.Vpl./rcv.Vo))+G*rcv.Vpl./(2*rcv.Vs); %stress 
Y0(3:rcv.dgf:end) = log(rcv.Vo./rcv.Vpl);      % log(theta Vo / L)
Y0(4:rcv.dgf:end) = log(rcv.Vpl*0.999./rcv.Vo); % log(V/Vo)

% initialize the function handle with
% set constitutive parameters
% yp=@(t,y) eqcycle_accelerationage_pinnedmegathrust2d_ode(t,y,rcv,K);
yp=@(t,y) eqcycle_accelerationslip_pinnedmegathrust2d_ode(t,y,rcv,K);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-7,'InitialStep',1e-3,'MaxStep',1e7); 
[ti,Yi]=ode45(yp,[0 tperturb(1)],Y0,options);
toc
t = ti;
Y = Yi;
% add a perturbation and continue with solution (problem is that right now
% stress perturbations have no effect, and need to implemented using
% velocity!
for i = 1:length(tperturb)
    tic
    Y0 = Y(end,:)';
    % v+ = v-*exp(deltau/a*sigma) -> transform to log(v/v0)
    Y0(4*top:rcv.dgf:bottom*rcv.dgf) = Y0(4*top:rcv.dgf:bottom*rcv.dgf) + deltau(i)./(rcv.a(top:bottom).*rcv.sigma(top:bottom));
    if i==length(tperturb)
        [ti,Yi]=ode45(yp,[0 tend-tperturb(i)],Y0,options);
    else
        [ti,Yi]=ode45(yp,[0 tperturb(i+1)-tperturb(i)],Y0,options);
    end
    toc
    
    t = [t;ti+tperturb(i)];
    Y = [Y;Yi];
end
% Velocities
V=repmat(rcv.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:rcv.dgf:rcv.N*rcv.dgf));
slip = Y(:,1:rcv.dgf:end);
% Maximum Velocity
Vmax=max(V,[],2);
return
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% make new colormap
% cspec=[cmap('steelblue',100,10,48);...
% %     cmap('skyblue',100,45,27);...
%     (cmap('skyblue',100,45,27));...
%     (cmap('lightgreen',10,59,32));...
%     flipud(cmap('orange',100,57,20));...
%     flipud(cmap('orangered',100,40,25))];
down = 50; %spatial downslampling

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(3,2,[1 3])
pcolor(t(1:down:end-1)./3.15e7,dipdist/1e3,log10(V(1:down:end,:)')), shading flat, hold on
plot([tperturb' tperturb']./3.15e7,get(gca,'YLim'),'k','LineWidth',2)
plot(get(gca,'XLim'),[dipdist(top) dipdist(top)]./1e3,'m-','Linewidth',1)
plot(get(gca,'XLim'),[dipdist(bottom) dipdist(bottom)]./1e3,'m-','Linewidth',1)
box on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
% colormap(cspec);
title(h,'log_{10} V'),ylabel('Depth (km)');
colormap(jet(32))
caxis([-11 -7])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,5);cla;
% plot(t(1:down:end-1)/3.15e7,(mean(V(1:down:end,top:bottom),2)./rcv.Vpl(1)),'m-','LineWidth',2), hold on
plot(t(1:down:end-1)/3.15e7,(max(V(1:down:end,top:bottom),[],2)./rcv.Vpl(1)),'m-','LineWidth',2), hold on
plot(t(1:down:end-1)/3.15e7,(max(V(1:down:end,topvw:botvw),[],2)./rcv.Vpl(1)),'b-','LineWidth',1)
plot([tperturb' tperturb']./3.15e7,get(gca,'YLim'),'k--','LineWidth',1)
% legend('VS_{perturb}','VW')
xlabel('Time (Yr)','Fontsize',15),ylabel('V/V_{pl}')
axis tight,grid on,ylim([0 2])
set(gca,'FontSize',15,'Color','none')

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(3,2,[2 4]);cla;
imagesc(1:down:length(t)-1,dipdist/1e3,log10(V(1:down:end,:)')), shading flat, hold on
plot([find(t==tperturb,1) find(t==tperturb,1)],get(gca,'Ylim'),'k','LineWidth',2)
plot(get(gca,'XLim'),[dipdist(top) dipdist(top)]./1e3,'m-','Linewidth',1)
plot(get(gca,'XLim'),[dipdist(bottom) dipdist(bottom)]./1e3,'m-','Linewidth',1)
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('Depth (km)');
caxis([-10 -8])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,6);cla;
plot(1:down:length(t)-1,(max(V(1:down:end,topvw:botvw),[],2)),'b-','LineWidth',1),hold on
plot(1:down:length(t)-1,(mean(V(1:down:end,top:bottom),2)),'m-','LineWidth',2)
xlabel('Time Steps','Fontsize',15),ylabel('V')
axis tight, grid on
plot([find(t==tperturb,1) find(t==tperturb,1)],get(gca,'Ylim'),'k--','LineWidth',1)
set(gca,'FontSize',15,'Color','none','YScale','log')


%% %% compute displacement field %%%%%%%%%%%%%%%%%%%%
ndisp = 20;
xplt = linspace(10e3,120e3,ndisp)';
[~,Gd]=rcv.displacementKernels([xplt,0.*xplt,0.*xplt],3);
[~,Gddeep]=rcvdeep.displacementKernels([xplt,0.*xplt,0.*xplt],3);
Gde = Gd(1:3:end,:);
Gdz = Gd(3:3:end,:);
discol = jet(ndisp);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%          Friction properties and model setup       %  %
figure(2);clf;
plot(rcv.a-rcv.b,dipdist./1e3,'k-','LineWidth',3)
hold on
plot(rcv.a,dipdist/1000,'b-','LineWidth',2)
plot(rcv.b,dipdist/1000,'r-','LineWidth',2)
legend('a-b','a','b')
set(legend,'box','off')
box on
axis tight
grid on
set(gca,'FontSize',15,'Color','none','YDir','reverse')
ylabel('Down-dip distance (km)')
xlabel('a-b')

%% cumulative slip plot

j = 1;
% initialize arrays for Yn = integrated array; 
% Vn = velocity extracted from array of derivatives
Yn = [];
Vn=[];
index =[];
Veq = 1e-3;

for i = 1:length(t)-1
    if i == 1
        index(j) = i;
        j = j+1;
    elseif i > 1 && max(V(i,:)) >=Veq  
        if ((t(i) - t(index(j-1))) > 3)
            index(j) = i;
            j = j+1;
        end
    else
        if ((t(i) - t(index(j-1))) > 40/365*3.15e7)
            index(j) = i;
            j = j+1;
        end
    end
end
Yn = Y(index',:);
Vn = V(index',:);


figure(3),clf
skip = 30;
count = 1;
for i = 1:(length(Yn(:,1))/1)
% for i = round(length(Yn(:,1))/2):length(Yn(:,1))
    if max(Vn(i,:)) >= Veq 
        plot(dipdist./1e3,Yn(i,1:rcv.dgf:end),'-','LineWidth',0.01,'Color',rgb('orangered')), hold on
    elseif (max(Vn(i,:)) <= 5*rcv.Vpl(1) && max(Vn(i,:)) >= rcv.Vpl(1))
        plot(dipdist./1e3,Yn(i,1:rcv.dgf:end),'-','LineWidth',1,'Color',rgb('forestgreen')), hold on
    else
        if count>=skip
            plot(dipdist./1e3,Yn(i,1:rcv.dgf:end),'-','LineWidth',1,'Color',rgb('steelblue')), hold on
            count = 1;
        else
            count = count+1;
        end
    end    
end
xlabel('Distance along fault (km)')
ylabel('Cumulative Slip (m)')
axis tight, grid off
set(gca,'FontSize',15,'Color','none')
