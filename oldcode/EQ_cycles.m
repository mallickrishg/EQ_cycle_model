% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                         %
%         E a r t h q u a k e   c y c l e s               %
%      o n   a   s t r i k e - s l i p   f a u l t        %
%                                                         %
% DESCRIPTION:                                            %
% Solves the governing equations for fault slip evolution % 
% using rate-and-state friction using the velocity in the %
% state vector and the acceleration in the state          %
% derivative.                                             %
%                                                         %
% Evaluates the evolution of slip using the boundary      %
% integral method.                                        %
%                                                         %
% AUTHOR:                                                 %
% Adapted from Sylvain and Valere's example               %
% by Rishav Mallick                                       %         
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear

minmax=@(x) [min(x(:)),max(x(:))];

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
Transition = 35e3;

% number of fault patches
ss.M = 300;

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
% stresses that are hydrostatic (2700-1000)kg/m3
% ss.sigmab = 1700*10*(ss.y3f+dz/2)*1e-6;
% ss.sigmab(ss.sigmab < 100) = 100;

ss.sigmab = 100.*ones(size(ss.y3f));
% ss.sigmab(1:ceil(5e3/dz)) = linspace(100,500,ceil(5e3/dz))';
% frictional parameters ( velocity-weakening friction, a-b < 0 )
ss.a = 1e-2*ones(size(ss.y3f));
ss.b = ss.a+5e-3*ones(size(ss.y3f));

% static friction coefficient
ss.mu0 = 0.6*ones(size(ss.y3f));

% characteristic weakening distance (m)
ss.L = 0.01*ones(size(ss.y3f));

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
top    = floor(5e3/(Transition/ss.M));
bottom = ceil(15e3/(Transition/ss.M));

% top1 = floor(5e3/(Transition/ss.M));
% bottom1 = ceil(15e3/(Transition/ss.M));
% top2 = floor(20e3/(Transition/ss.M));
% bottom2 = ceil(30e3/(Transition/ss.M));

% b_slope=(ss.b(top+1)-ss.a(1))/top;
% ss.b(1:top)      = b_slope*(1:top)+.7e-3;
ss.b(1:top)      = ss.a(1:top)     -1e-3*ones(top,1);
ss.b(bottom:end) = ss.a(bottom:end)-5e-3*ones(length(ss.a(bottom:end)),1);

% ss.b(1:top1)      = -ss.a(1:top1);
% ss.b(bottom1:top2) = ss.a(bottom1:top2)     -2e-3*ones(top2+1-bottom1,1);
% ss.b(bottom2:end) = -ss.a(bottom2:end);

% fitb = polyfit(ss.y3f,ss.b,8);
% ss.b = polyval(fitb,ss.y3f);
% ss.b(bottom:ceil(25e3/(Transition/ss.M))) = ss.b(bottom).*ones(size(ss.b(bottom:ceil(25e3/(Transition/ss.M))))) + ...
%     (ss.b(end)-ss.b(bottom-1))/top.*(bottom:ceil(25e3/(Transition/ss.M)))' + 2.11e-3;
% Fault Strength
ss.strength = ss.sigmab.*(ss.mu0+(ss.a-ss.b).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         N U M E R I C A L   S O L U T I O N          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% state parameters
ss.dgfF=4;

%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF,1);

% Fault patches
Y0(1:ss.dgfF:end) = zeros(size(ss.y3f));         % slip
Y0(2:ss.dgfF:end) = ss.strength;                 % stress
Y0(3:ss.dgfF:end) = log(ss.Vo./ss.V_plate);      % log(theta Vo / L)
Y0(4:ss.dgfF:end) = log(ss.V_plate*0.98./ss.Vo); % log(V/Vo)

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) earthquake_cycle_acceleration_ode(t,y,ss);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6); 
[t,Y]=ode45(yp,[0 1.5e10],Y0,options);
toc

% Velocities
V=repmat(ss.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:ss.dgfF:ss.M*ss.dgfF));

% Maximum Velocity
Vmax=max(V,[],2);

% Velocity in shallow creeping bit
%Vcreep = max(V(:,1:top-1),[],2);%/max(ss.V_plate);
% Vcreep = mean(V(:,1:top-1),2);

%% save data
save('results/result_data_test','-v7.3')
clear
