function [sol,exitflag,solver,par] = optimize(sys,d_val,options)

% Computes the steady-state optimum at a given set of disturbance d_val
% Written by Dinesh Krishnamoorthy, Jul 2019, NTNU
% Uses CasADi v3.5.1

import casadi.*

if nargin<3
    % Default options
    options.opts = struct('warn_initial_bounds',false, ...
        'print_time',false, ...
        'ipopt',struct('print_level',1)...
        );
    
    options.lbx = 0.*ones(size(sys.x));
    options.ubx = 1e10.*ones(size(sys.x));
    options.x0 = 1e-2.*ones(size(sys.x));
    
    options.lbu = 0.*ones(size(sys.u));
    options.ubu = 1e10.*ones(size(sys.u));
    options.u0 = 1e-2.*ones(size(sys.u));
end

lbx = options.lbx;
ubx = options.ubx;
x0 = options.x0;

lbu = options.lbu;
ubu = options.ubu;
u0 = options.u0;

assert(numel(sys.d)==numel(d_val))

w = {};
w0 = [];
lbw = [];
ubw = [];

g = {};
lbg = [];
ubg = [];

w = {w{:},sys.x,sys.u};
lbw = [lbw;lbx;lbu];
ubw = [ubw;ubx;ubu];
w0 = [w0;x0;u0];

J = sys.L; % Economic objective


g = {g{:},vertcat(sys.f)};
lbg = [lbg;zeros(numel(sys.f),1)];
ubg = [ubg;zeros(numel(sys.f),1)];

if ~isempty(sys.nlcon)
    assert(numel(sys.nlcon)==numel(sys.lb))
    assert(numel(sys.nlcon)==numel(sys.ub))
    
    g = {g{:},sys.nlcon};
    lbg = [lbg;sys.lb];
    ubg = [ubg;sys.ub];
end

nlp = struct('x',vertcat(w{:}),'p',vertcat(sys.d),'f',J,'g',vertcat(g{:}));
solver = nlpsol('solver','ipopt',nlp,options.opts);

sol = solver('x0',w0,'p',vertcat(d_val),'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
wf = full(sol.x);
xf = wf(1:numel(sys.x));

flag = solver.stats();
exitflag =  flag.success;

par.w0 = w0;
par.lbw = lbw;
par.ubw = ubw;
par.lbg = lbg;
par.ubg = ubg;
par.nlp = nlp;



