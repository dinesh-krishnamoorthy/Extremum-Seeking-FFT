function P_tot = windFarm6(u,Vinf)
% Wind turbine model from Marden et al. 2013, IEEE Trans CST

if nargin <2
    Vinf = 8;
end

par.D   = 80; % m
par.k   = 0.075; % blade roughness co-efficient. 0.04 for offshore
par.rho = 1.225; % kg/m3
par.x   = [0;400;800;0;400;800]; % location of the turbines
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

% Turbine 4
V(4) = par.Vinf;
Cp(4) = 4*u(4)*(1-u(4))^2;
P(4) = 0.5*par.rho*A*Cp(4)*V(4)^3;

% for turbine 5
c(4,5) = (par.D/(par.D + 2*par.k*(par.x(5) - par.x(4))))^2;
V(5) = par.Vinf*(1-2*sqrt((u(4)*c(4,5))^2));
Cp(5) = 4*u(5)*(1-u(5))^2;
P(5) = 0.5*par.rho*A*Cp(5)*V(5)^3;

% for turbine 6
c(4,6) = (par.D/(par.D + 2*par.k*(par.x(6) - par.x(4))))^2;
c(5,6) = (par.D/(par.D + 2*par.k*(par.x(6) - par.x(5))))^2;
V(6) = par.Vinf*(1-2*sqrt((u(4)*c(4,6) + u(5)*c(5,6))^2));
Cp(6) = 4*u(6)*(1-u(6))^2;
P(6) = 0.5*par.rho*A*Cp(6)*V(6)^3;

P_tot = sum(P);


