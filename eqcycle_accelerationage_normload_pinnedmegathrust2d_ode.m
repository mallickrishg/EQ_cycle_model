function [Yp]= eqcycle_accelerationage_normload_pinnedmegathrust2d_ode(~,Y,rcv,K,dsigma_dt)            
% function describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%        |       tau         |         
%    y = | log(theta Vo / L) |           
%        |   log( V / Vo )   |
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
damping = rcv.earthModel.G./rcv.Vs/2;
% Shear stress on faults
tauF = Y(2:rcv.dgf:end);

% normal stress
sigma = Y(5:rcv.dgf:end);

% friction coefiicient
% State variables
th   = Y(3:rcv.dgf:end);
th(rcv.pinnedPosition) = 0;

% Slip velocities 
V    = rcv.Vo.*exp(Y(4:rcv.dgf:end));
V(rcv.pinnedPosition) = 0;

f = rcv.mu0 + rcv.a.*Y(4:rcv.dgf:end) + rcv.b.*th;

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = V;

% Rate of state (rate of log(theta/theta0))
dth=(rcv.Vo.*exp(-th)-V)./rcv.l;
Yp(3:rcv.dgf:rcv.N*rcv.dgf) = dth;

% Acceleration (rate of log(V/Vo))
kv = K*(V - rcv.Vpl);
Yp(4:rcv.dgf:rcv.N*rcv.dgf) = (kv - rcv.b.*sigma.*dth - f.*dsigma_dt)./(rcv.a.*sigma+damping.*V);

% Shear stress rate on fault due to fault
Yp(2:rcv.dgf:rcv.N*rcv.dgf) = kv - damping.*V.*Yp(4:rcv.dgf:rcv.N*rcv.dgf);

% normal stress change rate
Yp(5:rcv.dgf:rcv.N*rcv.dgf) = dsigma_dt;
end

