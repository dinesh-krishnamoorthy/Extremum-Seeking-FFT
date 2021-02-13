function u_opt_true = true_plant_optimum(plant,d_val)
import casadi.*

% Default options
options.opts = struct('warn_initial_bounds',false, ...
    'print_time',false, ...
    'ipopt',struct('print_level',1)...
    );

options.lbx = 0.*ones(size(plant.x));
options.ubx = 1e10.*ones(size(plant.x));
options.x0 = [108.4165
    110.0000
    108.4165
    108.4165
    110.0000
    108.4165
    110.0000
    108.4165
    114.5819
    87.4620
    116.1290
    91.3666
    87.4620
    116.1290
    91.3666
    116.1290
    91.3666
    0.2];

options.lbu = 0.*ones(size(plant.u));
options.ubu = 1e10.*ones(size(plant.u));
options.u0 = 0.1.*ones(size(plant.u));

[sol,exitflag] = optimize(plant,d_val,options);
w_opt = full(sol.x);
u_opt_true = w_opt(numel(plant.x)+1:numel(plant.x)+numel(plant.u));

end
