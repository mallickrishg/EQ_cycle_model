% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                         %
%         E a r t h q u a k e   c y c l e s               %
%      o n   a   s t r i k e - s l i p   f a u l t        %
%                                                         %
% DESCRIPTION:                                            %
% Solves the governing equations for fault slip evolution % 
% on a purely velocity-strengthening fault using          %
% rate-and-state friction using the slip in the           %
% state vector and velocity in the state                  %
% derivative.                                             %
%                                                         %
% Evaluates the evolution of slip using the boundary      %
% integral method.                                        %
%                                                         %
% AUTHOR:                                                 %
% Rishav Mallick                                          %         
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear

% Rigidity (MPa)
G = 30e3;

% displacement kernel
d1=@(x2,y3,W) -1/pi.*(atan2(x2,y3+W/2) - atan2(x2,y3-W/2));

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

% Fault Meshes
y3 = 0e3; % Fault starting depth (m)
y2 = 0e3; % Fault horizontal position (m)

% Brittle-Ductile Tranisition Depth (m)
Transition = 20e3;

% number of fault patches (use the criteria dz < Lb =
% G*ss.L/ss.b/ss.sigmab) A good rule is to use dz < L/20
ss.M = 800;

dz     = Transition/ss.M;
fpoles = y3+(0:ss.M)'*dz;

% top of fault patches
ss.y3f = fpoles(1:end-1); 
% width of fault patches
Wf     = ones(ss.M,1)*dz; 

%% Create stress kernels for fault interactions

ss.K=zeros(ss.M,ss.M);   % Fault self stress

% Evaluate the stress at the center of the fault
for k=1:ss.M
    % Stress on faults from fault slip
    ss.K(:,k)=s12h(y2,ss.y3f+dz/2,y2,ss.y3f(k),Wf(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% effective confining pressure on fault (MPa) using depth dependent normal
ss.sigmab = 50.*ones(size(ss.y3f));

% frictional parameters 
ss.a = 1e-2*ones(size(ss.y3f));
% if you want to set the entire thing as VS
ss.b = ss.a - .2e-3*ones(size(ss.y3f));

% static friction coefficient
ss.mu0 = 0.6*ones(size(ss.y3f));

% characteristic weakening distance (m)
ss.L = 0.005*ones(size(ss.y3f));

% plate velocity (m/s)
ss.V_plate = 1e-9*ones(size(ss.y3f));

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(size(ss.y3f));

% shear wave speed (m/s)
ss.Vs = 3e3*ones(size(ss.y3f));

% radiation damping term
ss.damping=G./ss.Vs/2;

% Velocity-strengthening at top and bottom ( a-b > 0 )
% 5-15km velocity-weakening
top    = floor(10e3/(Transition/ss.M));
topup    = floor(7e3/(Transition/ss.M));
bottom = ceil(15e3/(Transition/ss.M));

ss.b(top:bottom)      = ss.a(top:bottom)     +8e-3;
% ss.b(bottom:end) = ss.a(bottom:end)-5e-3*ones(length(ss.a(bottom:end)),1);

% Fault Strength
ss.strength = ss.sigmab.*(ss.mu0+(ss.a-ss.b).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);
% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f\n',dz,G*ss.L(1)/ss.b(1)/ss.sigmab(1)) 
fprintf(1,'a/b = %.2f\n',ss.a(1)/ss.b(1))

% state parameters
ss.dgfF=4;
% add time of perturbation
tperturb = .245e10;
% stress perturbation
deltau = -.8;
% total duration of simulation
tend = .8e10;
%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF,1);

% Fault patches
Y0(1:ss.dgfF:end) = zeros(size(ss.y3f));         % slip
Y0(2:ss.dgfF:end) = ss.strength;                 % stress
Y0(3:ss.dgfF:end) = log(ss.Vo./ss.V_plate);      % log(theta Vo / L)
Y0(4:ss.dgfF:end) = log(ss.V_plate*0.99./ss.Vo); % log(V/Vo)

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) earthquake_cycle_acceleration_ode(t,y,ss);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-9,'InitialStep',1e-3,'MaxStep',3e6); 
[t1,Y1]=ode45(yp,[0 tperturb],Y0,options);
toc

% add a perturbation and continue with solution (problem is that right now
% stress perturbations have no effect, and need to implemented using
% velocity!
tic
Y0 = Y1(end,:)';
% v+ = v-*exp(deltau/a*sigma) -> transform to log(v/v0)
n = 80;
Y0(4*n:ss.dgfF:topup*ss.dgfF) = Y0(4*n:ss.dgfF:topup*ss.dgfF) + deltau./(ss.a(n:topup).*ss.sigmab(n:topup));

% Y0(4*bottom:ss.dgfF:end-n*ss.dgfF) = Y0(4*bottom:ss.dgfF:end-n*ss.dgfF) + deltau./(ss.a(bottom:end-n).*ss.sigmab(bottom:end-n));
% Y0(4:ss.dgfF:end) = Y0(4:ss.dgfF:end) + deltau./(ss.a.*ss.sigmab);

[t2,Y2]=ode45(yp,[0 tend-tperturb],Y0,options);
toc

t = [t1;t2+tperturb];
Y = [Y1;Y2];
% Velocities
V=repmat(ss.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:ss.dgfF:ss.M*ss.dgfF));

% Maximum Velocity
Vmax=max(V,[],2);


%% save data
save('results/result_shallow_transients','-v7.3')
clear
