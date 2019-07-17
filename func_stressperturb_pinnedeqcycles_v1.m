function [Y,V,t] = func_stressperturb_pinnedeqcycles_v1(patchfname,vwparams,perturbparams,tend)
% INPUTS
% patchfname - rcv path
% vwparams - structure that contains [a,b,top,bot] (does not use a - this is fixed at 0.01)
%          top - depth in m
%          bot - depth in m 
%                region from top:bot is pinned at does not take part in
%                earthquake cycle
% perturbparams - a,b,top,bottom,tperturb,deltau
%          top - depth in m
%          bot - depth in m 
%          tperturb - time in (s) typically 0.3e10 (give enough time for spin up (can be a vector)
%          deltau - stress change in MPa (use -0.5 for simuelue style transients) (must be a vector same size as tperturb)
%
% OUTPUTS - solved using rate-and-state friction and the ageing law
% Y - full state vector [slip, stress, log(theta Vo / L), log(V/Vo)]
% V - velocity vector in (m/s)
% t - time vector (in sec)
% Rishav Mallick, EOS, 2019

% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu
import unicycle.*
G = 30e3; nu = 0.25;
earthModel=unicycle.greens.okada92(G,nu);

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

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% create velocity weakening region
VW = zeros(rcv.N,1);
topvw = floor(vwparams.top/(rcv.W(1)));
botvw = ceil(vwparams.bot/(rcv.W(1)));
VW(topvw:botvw) = 1;
VW = logical(VW);
% effective confining pressure on fault (MPa) using depth dependent normal
rcv.sigma = 30.*ones(rcv.N,1);
% frictional parameters
rcv.a = vwparams.a.*ones(rcv.N,1);
rcv.b = vwparams.bvs.*ones(rcv.N,1); % for VS region
rcv.pinnedPosition = VW;

% static friction coefficient
rcv.mu0 = 0.6*ones(rcv.N,1);
% characteristic weakening distance (m)
rcv.l = 0.005.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 1e-9*ones(rcv.N,1);
% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);
% location of perturbation
top    = floor(perturbparams.top/(rcv.W(1)));
bottom = ceil(perturbparams.bot/(rcv.W(1)));
rcv.b(top:bottom) = perturbparams.b;
% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),G*rcv.l(1)/rcv.b(top+1)/rcv.sigma(1))
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(1)/rcv.b(1))

% numer of state parameters
rcv.dgf=4;
% add time of perturbation
tperturb = perturbparams.tperturb;
% stress perturbation
deltau = perturbparams.deltau;
% check if tperturb and deltau are the same length
if length(tperturb) ~= length(deltau)
    error('Stress Perturbations and Event Time vector need to be same length')
end
% total duration of simulation (taken as input)
% tend = .8e10;
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
options=odeset('Refine',1,'RelTol',3e-9,'InitialStep',1e-3,'MaxStep',3e6);
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
V = repmat(rcv.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:rcv.dgf:rcv.N*rcv.dgf));

end