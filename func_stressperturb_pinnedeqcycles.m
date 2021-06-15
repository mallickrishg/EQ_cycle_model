function [Y,V,t] = func_stressperturb_pinnedeqcycles(rcv,tperturb,deltau,top,bottom,tend,thetalaw)
% INPUTS
% rcv - unicycle fault
% tperturb - can be scalar or vector of time at which stress perturbations occur
% deltau - scalar/vector (same size as tperturb) of stress perturbations
% top,bottom - region to apply stress perturbation (provide as index for rcv)
% tend - total duration of simulation
% thetalaw (state evolution) - (1) ageing law ; (2) slip law
% OUTPUTS - solved using rate-and-state friction and the ageing law
% Y - full state vector [slip, stress, log(theta Vo / L), log(V/Vo)]
% V - velocity vector in (m/s)
% t - time vector (in sec)
% Rishav Mallick, EOS, 2019

% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu

G = 30e3; nu = 0.25;
earthModel=unicycle.greens.okada92(G,nu);


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

% minimum grid size
% fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),G*rcv.l(1)/max(rcv.b)/rcv.sigma(1))
% fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(1)/rcv.b(1))


% check if tperturb and deltau are the same length
if length(tperturb) ~= length(deltau)
    error('Stress Perturbations and Event Time vector need to be same length')
end

%% Initialize State Vector
Y0=zeros(rcv.N*rcv.dgf,1);
% For fault patches
Y0(1:rcv.dgf:end) = zeros(rcv.N,1);         % slip
Y0(2:rcv.dgf:end) = rcv.sigma.*(rcv.mu0+(rcv.a-rcv.b).*log(rcv.Vpl./rcv.Vo))+G*rcv.Vpl./(2*rcv.Vs); %stress
Y0(3:rcv.dgf:end) = log(rcv.Vo./rcv.Vpl);      % log(theta Vo / L)
Y0(4:rcv.dgf:end) = log(rcv.Vpl*0.999./rcv.Vo); % log(V/Vo)

% initialize the function handle with
% set constitutive parameters
if thetalaw == 1
    disp('Using Ageing Law')
    yp=@(t,y) eqcycle_accelerationage_pinnedmegathrust2d_ode(t,y,rcv,K);
else
    disp('Using Slip Law')
    yp=@(t,y) eqcycle_accelerationslip_pinnedmegathrust2d_ode(t,y,rcv,K);
end
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