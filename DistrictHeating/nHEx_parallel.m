function sys = nHEx_parallel(nHEx)

% Written by: Dinesh Krishnamoorthy, July 2020
% Uses CasADi v3.5.1 www.casadi.org

import casadi.*

% states
Ti = MX.sym('Ti',nHEx);   % outlet temperature of each cold split stream
T = MX.sym('T');          % total output temperature
The = MX.sym('The',nHEx); % Outlet temperatue of hot stream  
ul = MX.sym('ul',1);      % fraction through the last branch

% input
u = MX.sym('u',nHEx-1); % split ratios

% disturbances
T0 = MX.sym('T0');          % inlet temperature of the cold stream
w0 = MX.sym('w0');          % inlet flow rate of cold stream
Th = MX.sym('Th',nHEx);     % inlet tempreature of each hot stream
wh = MX.sym('wh',nHEx);     % inlet flow rate of each hot stream
UA = MX.sym('UA',nHEx);     % UA of each HEx

f1 = -Ti + T0 + (Th-T0)./(1/2+[u;ul].*w0.*(1./(2*wh)+1./UA));
f2 = -T + [u;ul]'*Ti;
f3 = -The + Th - [u;ul].*w0./wh.*(Ti - T0);
f4 = -ul + 1 - sum(u); % Mass balance

sys.x = vertcat(Ti,T,The,ul);              % states
sys.u = vertcat(u);                              % inputs
sys.d = vertcat(T0,w0,Th,wh,UA);           % disturbances
sys.y = sys.x;%vertcat(T0,Ti,Th,The,T);        % measurements
sys.f = vertcat(f1,f2,f3,f4);              % model equations
sys.L = -T;                             % objective function
sys.nlcon = [];
sys.lb = []; sys.ub = [];


