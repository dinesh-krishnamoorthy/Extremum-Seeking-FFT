function [xf,exitflag] = solvef(sys,d_val,u_in,options)

% Rootfinder function that computes the steady-state at a given input u_in
% Written by Dinesh Krishnamoorthy, Jul 2019, NTNU
% Uses CasADi v3.5.1

import casadi.*

if nargin<4
    % Default options
    options.opts = struct('warn_initial_bounds',false, ...
        'print_time',false, ...
        'ipopt',struct('print_level',1)...
        );
    
    options.lbx = 0.*ones(size(sys.x));
    options.ubx = 1e10.*ones(size(sys.x));
    options.x0 = 1e-2.*ones(size(sys.x));
end

lbx = options.lbx;
ubx = options.ubx;
x0 = options.x0;

assert(numel(sys.d)==numel(d_val))
assert(numel(sys.u)==numel(u_in))

w = {sys.x};

g = {vertcat(sys.f)};
lbg = zeros(numel(sys.f),1);
ubg = zeros(numel(sys.f),1);

nlp = struct('x',vertcat(w{:}),'p',vertcat(sys.u,sys.d),'f',0,'g',vertcat(g{:}));
solver = nlpsol('solver','ipopt',nlp,options.opts);

sol = solver('x0',x0,'p',vertcat(u_in,d_val),'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
wf = full(sol.x);
xf = wf(1:numel(sys.x));

flag = solver.stats();
exitflag =  flag.return_status;


