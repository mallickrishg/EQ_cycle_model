 function [Yp]= eqcycle_accelerationage_deepdriver2d_ode(~,Y,rcv,K,Kn,tau_inf,sigma_inf)            
% function odeViscoelastic describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%        |       tau         |         
%    y = | log(theta Vo / L) |           
%        |   log( V / Vo )   |
%        \       sigma         /
%
%

damping = rcv.earthModel.G./rcv.Vs/2;
% Shear stress on faults
%tauF = Y(2:rcv.dgf:end);

% State variables
th   = Y(3:rcv.dgf:end);

% confining stress
% tolsigma = 40;
sigma = Y(5:rcv.dgf:end);
% sigma(sigma<(rcv.sigma-tolsigma)) = (rcv.sigma(sigma<(rcv.sigma-tolsigma)) - tolsigma);
% sigma(sigma>(rcv.sigma+tolsigma)) = (rcv.sigma(sigma>(rcv.sigma-tolsigma)) + tolsigma);

% Slip velocities 
V    = rcv.Vo.*exp(Y(4:rcv.dgf:end));

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = V;

% Rate of state (rate of log(theta/theta0))
dth=(rcv.Vo.*exp(-th)-V)./rcv.l;
Yp(3:rcv.dgf:rcv.N*rcv.dgf) = dth;

% normal stress change
dsigdt = (sigma_inf + Kn*V);
Yp(5:rcv.dgf:rcv.N*rcv.dgf) = dsigdt;

% Acceleration (rate of log(V/Vo))
kv = tau_inf + K*V;
Yp(4:rcv.dgf:rcv.N*rcv.dgf) = (kv - rcv.b.*sigma.*dth - (rcv.mu0 + rcv.a.*Y(4:rcv.dgf:end) + rcv.b.*th).*dsigdt)./(rcv.a.*sigma + damping.*V);

% Shear stress rate on fault due to fault
Yp(2:rcv.dgf:rcv.N*rcv.dgf) = kv - damping.*V.*Yp(4:rcv.dgf:rcv.N*rcv.dgf);

end

