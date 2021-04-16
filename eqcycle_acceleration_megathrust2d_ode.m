 function [Yp]= eqcycle_acceleration_megathrust2d_ode(~,Y,rcv,K)            
% function odeViscoelastic describes the evolution of the ordinary
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

% State variables
th   = Y(3:rcv.dgf:end);

% Slip velocities 
V    = rcv.Vo.*exp(Y(4:rcv.dgf:end));

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = V;

% Rate of state (rate of log(theta/theta0))
dth=(rcv.Vo.*exp(-th)-V)./rcv.l;
Yp(3:rcv.dgf:rcv.N*rcv.dgf) = dth;

% Acceleration (rate of log(V/Vo))
kv = K*(V - rcv.Vpl);
Yp(4:rcv.dgf:rcv.N*rcv.dgf) = (kv - rcv.b.*rcv.sigma.*dth)./(rcv.a.*rcv.sigma+damping.*V);

% Shear stress rate on fault due to fault
Yp(2:rcv.dgf:rcv.N*rcv.dgf) = kv - damping.*V.*Yp(4:rcv.dgf:rcv.N*rcv.dgf);

end

