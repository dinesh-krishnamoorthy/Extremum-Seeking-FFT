function P = windFarm(u,Vinf)
% Wind turbine model from Marden et al. 2013, IEEE Trans CST

if nargin <2
    Vinf = 8;
end
n = 3;

par.D   = 80; % m
par.k   = 0.075; % blade roughness co-efficient. 0.04 for offshore
par.rho = 1.225; % kg/m3
par.x   = [0;400;800]; % location of the turbines
par.Vinf= Vinf; % wind speed in m/s

A = pi*par.D^2/4;
% Turbine 1

V(1) = par.Vinf;
Cp(1) = 4*u(1)*(1-u(1))^2;
P(1) = 0.5*par.rho*A*Cp(1)*V(1)^3;

% for turbine 2
c(1,2) = (par.D/(par.D + 2*par.k*(par.x(2) - par.x(1))))^2;
V(2) = par.Vinf*(1-2*sqrt((u(1)*c(1,2))^2));
Cp(2) = 4*u(2)*(1-u(2))^2;
P(2) = 0.5*par.rho*A*Cp(2)*V(2)^3;

% for turbine 3
c(1,3) = (par.D/(par.D + 2*par.k*(par.x(3) - par.x(1))))^2;
c(2,3) = (par.D/(par.D + 2*par.k*(par.x(3) - par.x(2))))^2;
V(3) = par.Vinf*(1-2*sqrt((u(1)*c(1,3) + u(2)*c(2,3))^2));
Cp(3) = 4*u(3)*(1-u(3))^2;
P(3) = 0.5*par.rho*A*Cp(3)*V(3)^3;

P = sum(P);


