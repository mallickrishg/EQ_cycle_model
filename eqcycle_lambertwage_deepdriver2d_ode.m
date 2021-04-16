 function [Yp]= eqcycle_lambertwage_deepdriver2d_ode(~,Y,rcv,K,tau_inf)            
% function odeViscoelastic describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%        |       tau         |         
%    y = | log(theta Vo / L) |           
%        \       ...         /
%
% Instead of integrating numerically the aging law
%
%    d theta / dt = 1 - V theta / L
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / L)
%
% Considering the following relationships
%
%    d phi = d theta / theta
%
%    d theta = d phi theta = d phi exp(phi) L / Vo
%
%    d theta / dt = d phi exp(phi) / dt L / Vo
%
%    d phi exp(phi) / dt L / Vo = 1 - V theta / L = 1 - V exp(phi) / Vo
%
% we obtain the evolution law for the new variable
% 
%    d phi / dt = ( Vo exp(-phi) - V ) / L
%


% Shear stress on faults
tauF = Y(2:rcv.dgf:end);

% State variables
th   = Y(3:rcv.dgf:end);

% Slip velocities 
V    = (2.*rcv.Vs.*rcv.a.*rcv.sigma./rcv.earthModel.G).*...
      Lambert_W(rcv.earthModel.G*rcv.Vo./(2*rcv.Vs.*rcv.a.*rcv.sigma).*...
      exp((tauF-rcv.mu0.*rcv.sigma-rcv.sigma.*rcv.b.*th)./(rcv.sigma.*rcv.a)));

% Rate of state (rate of log(theta/theta0))
dth=(rcv.Vo.*exp(-th)-V)./rcv.l;

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = V;
% Shear stress rate on fault due to fault and shear zones
Yp(2:rcv.dgf:rcv.N*rcv.dgf) = K*V + tau_inf;
% Rate of state
Yp(3:rcv.dgf:rcv.N*rcv.dgf) = dth;

end

