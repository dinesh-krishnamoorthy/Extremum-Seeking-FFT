function sys = nHEx_parallel(nHEx)

% Written by: Dinesh Krishnamoorthy, July 2020
% Uses CasADi v3.5.1 www.casadi.org

import casadi.*

% states
Tc_out = MX.sym('Tc_out',nHEx);   % outlet temperature of each cold split stream
T = MX.sym('T');          % total output temperature
Th_out = MX.sym('Th_out',nHEx); % Outlet temperatue of hot stream
ul = MX.sym('ul',1);      % fraction through the last branch

% input
u = MX.sym('u',nHEx-1); % split ratios

% disturbances
Tc_in = MX.sym('Tc_in');          % inlet temperature of the cold stream
wc = MX.sym('wc');          % inlet flow rate of cold stream
Th_in = MX.sym('Th',nHEx);     % inlet tempreature of each hot stream
wh = MX.sym('wh',nHEx);     % inlet flow rate of each hot stream
UA = MX.sym('UA',nHEx);     % UA of each HEx

cp = 4.2;

DT2 = Th_in - Tc_out;
DT1 = Th_out - Tc_in;
DT_lm = (DT1.*DT2.*(DT1+DT2)/2).^(1/3); % Chen's approximation of the log mean temeprature

% DT_lm = (Th_in+Th_out)./2 - (Tc_in+Tc_out)./2; % Arithmetic mean

f1 = -Tc_out + Tc_in + UA.*DT_lm./([u;ul].*wc.*cp);
f2 = -T + [u;ul]'*Tc_out;
f3 = -Th_out + Th_in - [u;ul].*wc./wh.*(Tc_out - Tc_in);
f4 = -ul + 1 - sum(u); % Mass balance

sys.x = vertcat(Tc_out,T,Th_out,ul);              % states
sys.u = vertcat(u);                              % inputs
sys.d = vertcat(Tc_in,wc,Th_in,wh,UA);           % disturbances
sys.y = sys.x;
sys.f = vertcat(f1,f2,f3,f4);              % model equations
sys.L = -T;                             % objective function
sys.nlcon = [];
sys.lb = []; sys.ub = [];


