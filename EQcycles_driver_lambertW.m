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

patchfname='flt_2D_testing_v2.seg';

patchdeepname='flt_2D_deep_v2.seg';
%% 1.generate the fault geometry
dip1 = 30;
dip2 = 8;
dip3 = 30;
v_plate = 1;
fault_width = 15*1e3; %meters - along-dip width
npatch = 200;
patch_width = fault_width/npatch;
% patchfname='flt_2D_testing.seg';

% fileID = fopen(patchfname,'w');
% fprintf(fileID,'%s\n','#patch file generated automatically - 2D ramp model, constant patch size');
% fprintf(fileID,'%s\n',        '# n  Vpl    x1      x2   x3   Length  Width   Strike  Dip  Rake      L0     W0    qL qW');
% fprintf(fileID,'%s %.9f %s %.9f %s %.9f %s %.9f %s\n', '  1  ', v_plate, ...
%     ' -250e5       0    0e3    500e5 ', fault_width/1, ' 0 ', dip1, ' 90   500e5  ',patch_width,'  1  1 ');
% fprintf(fileID,'%s %.9f %s %.9f %.9f %s %.9f %s %.9f %s %.9f %s\n', '  1  ', v_plate, ...
%     ' -250e5    ' ,  fault_width*cosd(dip1)  ,  fault_width*sind(dip1)  ,'  500e5 ', fault_width/1, ' 0 ', dip2, ' 90   500e5  ',patch_width,'  1  1 ');
% fprintf(fileID,'%s %.9f %s %.9f %.9f %s %.9f %s %.9f %s %.9f %s\n', '  1  ', v_plate, ...
%     ' -250e5    ' ,  fault_width*(cosd(dip1)+cosd(dip2))  ,  fault_width*(sind(dip1)+sind(dip2))  ,'  500e5 ', fault_width/1, ' 0 ', dip3, ' 90   500e5  ',patch_width,'  1  1 ');
% 
% fclose(fileID);


% patchdeepname='flt_2D_deep.seg';
% fileID = fopen(patchdeepname,'w');
% fprintf(fileID,'%s\n','#patch file generated automatically - 2D ramp model, constant patch size');
% fprintf(fileID,'%s\n',        '# n  Vpl    x1      x2   x3   Length  Width   Strike  Dip  Rake      L0     W0    qL qW');
% fprintf(fileID,'%s %.9f %s %.9f %.9f %s %.9f %s %.9f %s %.9f %s\n', '  1  ', v_plate, ...
%     ' -250e5   ',fault_width*(cosd(dip1)+cosd(dip2)+cosd(dip3)),fault_width*(sind(dip1)+sind(dip2)+sind(dip3)),'  500e5 ', 1e8, ' 0 ', dip2, ' 90   500e5  ',1e8,'  1  1 ');
% fclose(fileID);
%% 2. compute greens functions and stress kernels
% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu
G = 30e3; nu = 0.25;
earthModel=greens.okada92(G,nu);
% rcv is the fault
rcv = unicycle.geometry.receiver(patchfname,earthModel);

T = 50e3;
[~,rcvdeep] = create_faults_exp(fault_width*(cosd(dip1)+cosd(dip2)+cosd(dip3)),fault_width*(sind(dip1)+sind(dip2)+sind(dip3)),0,2,T,'dummy.seg',patchdeepname);
% rcvdeep = geometry.receiver(patchdeepname,earthModel);
[~,~,~,Kdeep,~,~]=earthModel.tractionKernels(rcvdeep,rcv);

% self-stress kernels on the fault
% load stress kernels if they aren't already in the workspace
if ~exist('K','var')
    if exist('stresskernels/Kv2.mat','file')
        disp('Loading Kv2 from available file')
        load stresskernels/Kv2.mat
        if length(K)~=rcv.N
            disp('reloading appropriate kernels and fault file')                        
            [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
            save('stresskernels/Kv2.mat','K')
        end
    else
        disp('No Kv2 in folder, will create traction kernels')
        [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
        save('stresskernels/Kv2.mat','K')
    end
end
% define dip-distance in (m)
dipdist = zeros(rcv.N,1);
for i = 1:rcv.N
    if i==1
        dipdist(i) = rcv.W(1)/2;
    else
        dipdist(i) = dipdist(i-1) + rcv.W(i);
    end
end
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
fault_width = max(dipdist);
% create velocity weakening region
VW = zeros(rcv.N,1);
topvw = floor(10e3/(fault_width/rcv.N));
botvw = ceil(30e3/(fault_width/rcv.N));
VW(topvw:botvw) = 1;
VW = logical(VW);
% make VW region pinned
% rcv.pinnedPosition = VW;
% effective confining pressure on fault (MPa) using depth dependent normal
rcv.sigma = 30.*ones(rcv.N,1);

% frictional parameters 
rcv.a = 1e-2.*ones(rcv.N,1);
% if you want to set the entire thing as VS
rcv.b = rcv.a - 5e-3;
rcv.b(VW) = rcv.a(VW) + 4e-3;
% static friction coefficient
rcv.mu0 = 0.1*ones(rcv.N,1);
% numer of state parameters
rcv.dgf=4;
% characteristic weakening distance (m)
rcv.l = 0.01.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 1e-9*ones(rcv.N,1);
tau_inf = Kdeep*(rcvdeep.Vpl.*rcv.Vpl(1));
% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);

% location and rheology of perturbation
top    = floor(0.1e3/(fault_width/rcv.N));
bottom = ceil(2e3/(fault_width/rcv.N));
rcv.b(top:bottom) = rcv.a(top:bottom) - 8e-3;


% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),G/(1-nu)*rcv.l(topvw+1)/rcv.b(topvw+1)/rcv.sigma(topvw+1)) 
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(top+1)/rcv.b(top+1))
fprintf(1,'a/b in VW regions = %.2f\n',rcv.a(topvw+1)/rcv.b(topvw+1))
% total duration of simulation
tend = 200*3.15e7;
%% Initialize State Vector
Y0=zeros(rcv.N*rcv.dgf,1);
% Fault patches
Y0(1:rcv.dgf:end) = zeros(rcv.N,1);         % slip
Y0(2:rcv.dgf:end) = rcv.sigma.*(rcv.mu0+(rcv.a-rcv.b).*log(rcv.Vpl./rcv.Vo))+G*rcv.Vpl./(2*rcv.Vs); %stress 
Y0(3:rcv.dgf:end) = log(rcv.Vo./rcv.Vpl);      % log(theta Vo / L)

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) eqcycle_lambertwage_deepdriver2d_ode(t,y,rcv,K,tau_inf);
% yp=@(t,y) eqcycle_acceleration_megathrust2d_ode(t,y,rcv,K);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-7,'InitialStep',1e-3,'MaxStep',1e7); 
[ti,Yi]=ode45(yp,[0 tend],Y0,options);
toc
t = ti;
Y = Yi;

% Compute the instantaneous derivative
Yp=zeros(length(t)-1,size(Y,2));
for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
end
V=Yp(:,1:rcv.dgf:end);
slip = Y(:,1:rcv.dgf:end);
% Maximum Velocity
Vmax=max(V,[],2);
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

down = 50; %spatial downslampling

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(3,2,[1 3])
pcolor(t(1:down:end-1)./3.15e7,dipdist/1e3,log10(V(1:down:end,:)')), shading flat, hold on
plot(get(gca,'XLim'),[dipdist(top) dipdist(top)]./1e3,'k-','Linewidth',1)
plot(get(gca,'XLim'),[dipdist(bottom) dipdist(bottom)]./1e3,'k-','Linewidth',1)
box on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
% colormap(cspec);
title(h,'log_{10} V'),ylabel('Depth (km)');
colormap(jet(40))
caxis([-13 0])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,5);cla;
plot(t(1:down:end-1)/3.15e7,(max(V(1:down:end,top:bottom),[],2)./rcv.Vpl(1)),'k-','LineWidth',3), hold on
plot(t(1:down:end-1)/3.15e7,(max(V(1:down:end,topvw:botvw),[],2)./rcv.Vpl(1)),'r-','LineWidth',2)
xlabel('Time (Yr)','Fontsize',15),ylabel('V/V_{pl}')
axis tight,grid on,
% ylim([0 10])
set(gca,'FontSize',15,'Color','none','YScale','log')
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(3,2,[2 4]);cla;
imagesc(1:down:length(t)-1,dipdist/1e3,log10(V(1:down:end,:)')), shading flat, hold on

plot(get(gca,'XLim'),[dipdist(top) dipdist(top)]./1e3,'k-','Linewidth',1)
plot(get(gca,'XLim'),[dipdist(bottom) dipdist(bottom)]./1e3,'k-','Linewidth',1)
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('Depth (km)');
caxis([-13 0])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,6);cla;
plot(1:down:length(t)-1,(mean(V(1:down:end,top:bottom),2)),'k-','LineWidth',2), hold on
plot(1:down:length(t)-1,(mean(V(1:down:end,topvw:botvw),2)),'r-','LineWidth',2)
xlabel('Time Steps','Fontsize',15),ylabel('V')
axis tight, grid on

set(gca,'FontSize',15,'Color','none','YScale','log')



%% % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Slip coloured by velocity               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
figure(4),clf;set(gcf,'name','Time Evolution')
down = 1;
downt = 10;
ydata=repmat(dipdist(1:down:end)./1e3,[1,size(t(1:downt:end-1))])';
xdata=slip(1:downt:end-1,1:down:end);
zdata=log10(V(1:downt:end,1:down:end));

xsample=linspace(min(min(xdata(:))),max(max(xdata(:))),1e3);

test1=griddata(xdata(:),ydata(:),zdata(:),xsample,dipdist(1:down:end)./1e3);

pcolor(xsample,dipdist(1:down:end)./1e3,test1), shading flat

set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
caxis([-12 0]);
colormap(jet(24));
title(h,'Log_{10} V (m/s)')
xlabel('Accumulated slip (m)')
ylabel('Down-dip distance (km)');
set(gca,'FontSize',15)