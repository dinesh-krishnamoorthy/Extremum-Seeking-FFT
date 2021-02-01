clear
clc

FileName = mfilename('fullpath');
[directory,~,~] = fileparts(FileName);
[parent,~,~] = fileparts(directory);
addpath([parent '/functions'])

nL = 8;

L = 128/nL;
Fs = 128;
T = 1;

f = [6;11;17;23;31;39;47];
u0 = [0.1;0.1;0.1;0.1;0.1;0.1;0.1];
Ki = 0.1.*[0.00005;0.00005;0.00005;0.00005;0.00005;...
    0.00005;0.00005];
for k = 1:4000
    
    % Add dither
    a = [0.003;0.002;0.002;0.002;0.003;0.003;0.002];
    s = a.*sin(2*pi.*f.*k.*T./Fs);
    
    u = u0(:,k) + s;
    
    % Simulate
    J(k) = simulate_plant(u,k);
    
    sim.s(:,k) = s;
    sim.u(:,k) = u;
    sim.u0(:,k) = u0(:,k);
    sim.J(k) = J(k);
    
    % -------------
    
    if k>nL*L
        yL = sim.J(:,k-nL*L:k-1);
        uL = sim.u(:,k-nL*L:k-1);
        
        [Ymag,Yph,freq] = FFT(yL,nL,L,Fs);
        [U1mag,U1ph] = FFT(uL(1,:),nL,L,Fs);
        [U2mag,U2ph] = FFT(uL(2,:),nL,L,Fs);
        [U3mag,U3ph] = FFT(uL(3,:),nL,L,Fs);
        [U4mag,U4ph] = FFT(uL(4,:),nL,L,Fs);
        [U5mag,U5ph] = FFT(uL(5,:),nL,L,Fs);
        [U6mag,U6ph] = FFT(uL(6,:),nL,L,Fs);
        [U7mag,U7ph] = FFT(uL(7,:),nL,L,Fs);
        
        i1 = find(freq==f(1));
        i2 = find(freq==f(2));
        i3 = find(freq==f(3));
        i4 = find(freq==f(4));
        i5 = find(freq==f(5));
        i6 = find(freq==f(6));
        i7 = find(freq==f(7));
       
        sim.yph(1,k) = Yph(i1);
        sim.yph(2,k) = Yph(i2);
        sim.yph(3,k) = Yph(i3);
        sim.yph(4,k) = Yph(i4);
        sim.yph(5,k) = Yph(i5);
        sim.yph(6,k) = Yph(i6);
        sim.yph(7,k) = Yph(i7);
        
        sim.U1ph(k) = U1ph(i1);
        sim.U2ph(k) = U2ph(i2);
        sim.U3ph(k) = U3ph(i3);
        sim.U4ph(k) = U4ph(i4);
        sim.U5ph(k) = U5ph(i5);
        sim.U6ph(k) = U6ph(i6);
        sim.U7ph(k) = U7ph(i7);
        
        sim.Ju(1,k) = 2/sign(sim.yph(1,k).*sim.U1ph(k))*Ymag(i1)/U1mag(i1);
        sim.Ju(2,k) = 2/sign(sim.yph(2,k).*sim.U2ph(k))*Ymag(i2)/U2mag(i2);
        sim.Ju(3,k) = 2/sign(sim.yph(3,k).*sim.U3ph(k))*Ymag(i3)/U3mag(i3);
        sim.Ju(4,k) = 2/sign(sim.yph(4,k).*sim.U4ph(k))*Ymag(i4)/U4mag(i4);
        sim.Ju(5,k) = 2/sign(sim.yph(5,k).*sim.U5ph(k))*Ymag(i5)/U5mag(i5);
        sim.Ju(6,k) = 2/sign(sim.yph(6,k).*sim.U6ph(k))*Ymag(i6)/U6mag(i6);
        sim.Ju(7,k) = 2/sign(sim.yph(7,k).*sim.U7ph(k))*Ymag(i7)/U7mag(i7);
        
        u0(:,k+1) = u0(:,k) + Ki.*sim.Ju(:,k);
    else
        u0(:,k+1) = u0(:,k);
        
    end
end


%% plotting

figure(16)
clf
subplot(211)
hold all
plot(sim.J,'k-','linewidth',2)
ylabel('$T$','interpreter','latex')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
box on
grid on

subplot(212)
hold all
plot(sim.u0','linewidth',2)
ylabel('$u$','interpreter','latex')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
box on
grid on

%%

function J_meas = simulate_plant(u,k)

nHEx = 8;
plant = nHEx_parallel(nHEx);

par.T0=60;
par.w0=100;
par.Th = [120;120;120;120;120;120;120;120];
par.wh = [30;50;45;30;50;45;50;45];
par.UA = [50;55;60;50;55;60;55;60];

if k<2000
    par.Th = [120;120;120;120;120;120;120;120];
else
    par.Th = [150;120;120;120;120;120;120;120];
end
    
d_val = vertcat(par.T0,par.w0,par.Th,par.wh,par.UA);


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

xf = solvef(plant,d_val,u,options);
J_meas = xf(nHEx+1); 

end


function u_opt_true = plant_optimum(plant,d_val)

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
    114.5819
    87.4620
    116.1290
    91.3666];

options.lbu = 0.*ones(size(plant.u));
options.ubu = 1e10.*ones(size(plant.u));
options.u0 = [0.2016
    0.5323
    0.2661];


Ti_max = [120;120;120];
plant.nlcon = vertcat(plant.nlcon,plant.x(1),plant.x(2),plant.x(3));
plant.lb = [plant.lb;0;0;0];
plant.ub = [plant.ub;Ti_max];



[sol,exitflag] = optimize(plant,d_val,options);
w_opt = full(sol.x);
u_opt_all = w_opt(numel(plant.x)+1:numel(plant.x)+numel(plant.u));
u_opt_true = u_opt_all(1:2);

end
