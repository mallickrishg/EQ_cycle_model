clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*



%% 2. compute greens functions and stress kernels
% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu
G = 30e3; nu = 0.25;
earthModel=greens.okada92(G,nu);

patchfname='flt_2D_testing.seg';
% rcv is the fault
rcv = geometry.receiver(patchfname,earthModel);
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
rcv.mu0 = 0.2*ones(rcv.N,1);
% numer of state parameters
rcv.dgf=4;
% characteristic weakening distance (m)
rcv.l = 0.01.*ones(rcv.N,1);
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

% add time of perturbation
tperturb = 100*3.15e7;
% stress perturbation
deltau = [-0.5];
% total duration of simulation
tend = 500*3.15e7;

%% Vary parameter of interest and SOLVE (will be saved to directory)
thetalaw = 1; % 1 - ageing law, 2 - slip law
down = 30;

% vary parameter of choice
nexp = 10;
lvec = logspace(-3,-1,nexp);

for i = 1:nexp
    disp(['Iteration number ' num2str(i)])
 
    rcv.l = lvec(i).*ones(rcv.N,1);
    
    [Y,V,t] = func_stressperturb_pinnedeqcycles(rcv,tperturb,deltau,top,bottom,tend,thetalaw);

    % downsample spatially V and save it
    Vdown = V(:,1:down:end);
    % save files
    odir = '~/Documents/results_stressperturb_eqcycles';
    flist = dir([odir '/*.mat']);
    fnum = length(flist)+1;
    save([odir '/stressperturb_eqcycles_' num2str(fnum)],'rcv','tperturb','deltau','Vdown','t','-v7.3')
end
%% transfer all .mat to a new folder
flist = dir([odir '/run*']);
fnum = length(flist);
fol_out = [odir '/run' num2str(fnum+1)];
mkdir(fol_out)
movefile([odir '/*.mat'],fol_out);
%% visualize results

odatadir = ['~/Documents/results_stressperturb_eqcycles/run' num2str(fnum+1)];
figodir = ['Figures/perturbations/run' num2str(fnum+1)];
mkdir(figodir)

for i = 1:nexp
    datafnum = i;
    load([odatadir '/stressperturb_eqcycles_' num2str(datafnum) '.mat'])
    
    figure(1);clf
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    plot_stressperturb_pinnedeqcycles(rcv,Vdown,down,t)
    print([figodir '/pinnedeqcycle_deltau_' num2str(i)],'-djpeg','-r100')
    
    figure(2),clf
    
    %%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %          Friction properties and model setup       %  %
    plot(rcv.a-rcv.b,diag(rcv.xc*rcv.dv')./1e3,'k-','LineWidth',3)
    hold on
    plot(rcv.a,diag(rcv.xc*rcv.dv')/1000,'b-','LineWidth',2)
    plot(rcv.b,diag(rcv.xc*rcv.dv')/1000,'r-','LineWidth',2)
    legend('a-b','a','b')
    set(legend,'box','off')
    box on
    axis tight
    grid on
    set(gca,'FontSize',12,'Color','none')
    ylabel('Down-dip distance (km)')
    xlabel('a-b')
    print([figodir '/modelsetup_deltau_' num2str(i)],'-djpeg','-r100')
end