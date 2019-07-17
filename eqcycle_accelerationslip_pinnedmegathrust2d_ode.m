function [Yp]= eqcycle_accelerationslip_pinnedmegathrust2d_ode(~,Y,rcv,K)            
% function describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%        |       tau         |         
%    y = | log(theta Vo / L) |           
%        |   log( V / Vo )   |
%        \       ...         /
%
% we use the slip law to describe state parameter evolution

damping = rcv.earthModel.G./rcv.Vs/2;
% Shear stress on faults
tauF = Y(2:rcv.dgf:end);

% State variables
th   = Y(3:rcv.dgf:end);
th(rcv.pinnedPosition) = 0;

% Slip velocities 
V    = rcv.Vo.*exp(Y(4:rcv.dgf:end));
V(rcv.pinnedPosition) = 0;
zeta = Y(4:rcv.dgf:end);
zeta(rcv.pinnedPosition) = -Inf;

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = V;

% Rate of state (rate of log(theta Vo/L))
% dth = -(V./rcv.l).*(log(V./rcv.Vo) + th);
dth = -(rcv.Vo./rcv.l).*(zeta.*exp(zeta) + exp(zeta).*th);
dthetapinned = (rcv.Vo./rcv.l);
dth(rcv.pinnedPosition) = dthetapinned(rcv.pinnedPosition);
Yp(3:rcv.dgf:rcv.N*rcv.dgf) = dth;

% Acceleration (rate of log(V/Vo))
kv = K*(V - rcv.Vpl);
Yp(4:rcv.dgf:rcv.N*rcv.dgf) = (kv - rcv.b.*rcv.sigma.*dth)./(rcv.a.*rcv.sigma+damping.*V);

% Shear stress rate on fault due to fault
Yp(2:rcv.dgf:rcv.N*rcv.dgf) = kv - damping.*V.*Yp(4:rcv.dgf:rcv.N*rcv.dgf);

end

