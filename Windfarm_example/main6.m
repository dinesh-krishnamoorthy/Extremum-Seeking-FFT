clear
clc

FileName = mfilename('fullpath');
[directory,~,~] = fileparts(FileName);
[parent,~,~] = fileparts(directory);
addpath([parent '/functions'])

nL = 4;
L = 128/nL;
Fs = 128;
T = 1;

f = [6;17;31;39;47;11];  % Excitation frequency for the 6 wind turbines
u0 = [0.33;0.33;0.33;0.33;0.33;0.33]-0.03; % Initial input for the 6 turbines
Ki = [0.0005;0.0005;0.0005;0.00035;0.0005;0.0005]; % integral gain

for k = 1:1200
    
    sim.t(k) = k*T;
    Vinf = 8;  % wind velocity
  
    % Add dither
    a = 0.003;                      % dither amplitude
    s = a.*sin(2*pi.*f.*k.*T./Fs);
    u = u0(:,k) + s;
    
    % Simulate Wind farm with 6 turbines
    J(k) = windFarm6(u,Vinf)*1e-6 + 2e-3*randn(1,1);
    
    sim.s(:,k) = s;
    sim.u(:,k) = u;
    sim.u0(:,k) = u0(:,k);
    sim.J(k) = J(k);
    
    if k>nL*L
        yL = sim.J(:,k-nL*L:k-1);
        uL = sim.u(:,k-nL*L:k-1);
        
        % FFT of the cost signal
        [Ymag,Yph,freq,y0] = FFT(yL,nL,L,Fs);
        
        % FFT of the 6 input signals
        [U1mag,U1ph,freq,u10] = FFT(uL(1,:),nL,L,Fs);
        [U2mag,U2ph,freq,u20] = FFT(uL(2,:),nL,L,Fs);
        [U3mag,U3ph,freq,u30] = FFT(uL(3,:),nL,L,Fs);
        [U4mag,U4ph,freq,u40] = FFT(uL(4,:),nL,L,Fs);
        [U5mag,U5ph,freq,u50] = FFT(uL(5,:),nL,L,Fs);
        [U6mag,U6ph,freq,u60] = FFT(uL(6,:),nL,L,Fs);
        
        sim.Ymag(:,k) = Ymag;
        
        % pick out the frequncies of interest
        for i = 1:6
            ind(i) = find(freq==f(i));
        end
        
        sim.yph(:,k) = Yph(ind);
        
        sim.U1ph(k) = U1ph(ind(1));
        sim.U2ph(k) = U2ph(ind(2));
        sim.U3ph(k) = U3ph(ind(3));
        sim.U4ph(k) = U4ph(ind(4));
        sim.U5ph(k) = U5ph(ind(5));
        sim.U6ph(k) = U6ph(ind(6));
        
        sim.Ju(1,k) = sign(sim.yph(1,k).*sim.U1ph(k))*Ymag(ind(1))/U1mag(ind(1));
        sim.Ju(2,k) = sign(sim.yph(2,k).*sim.U2ph(k))*Ymag(ind(2))/U2mag(ind(2));
        sim.Ju(3,k) = sign(sim.yph(3,k).*sim.U3ph(k))*Ymag(ind(3))/U3mag(ind(3));
        sim.Ju(4,k) = sign(sim.yph(4,k).*sim.U4ph(k))*Ymag(ind(4))/U4mag(ind(4));
        sim.Ju(5,k) = sign(sim.yph(5,k).*sim.U5ph(k))*Ymag(ind(5))/U5mag(ind(5));
        sim.Ju(6,k) = sign(sim.yph(6,k).*sim.U6ph(k))*Ymag(ind(6))/U6mag(ind(6));
        
        u0(:,k+1) = u0(:,k) + 2*Ki.*sim.Ju(:,k);
    else
        u0(:,k+1) = u0(:,k);
        
    end
end

%%
figure(123)
time = nL*L:k;
mesh(time,freq/Fs,sim.Ymag(:,time),'linewidth',1,'EdgeColor','interp')
ylabel('Frequency','interpreter','latex')
xlabel('Time','interpreter','latex')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
c = colorbar('southoutside');
c.Label.String = '$2|\mathcal{J}(\omega)|$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';


figure(16)
clf
subplot(211)
hold all
plot(sim.J,'k-','linewidth',2)
ylabel('$P $ [MW]','interpreter','latex')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
box on
grid on
ylim([3.45,3.7])

subplot(212)
hold all
plot(sim.u0','linewidth',2)
ylabel('$u$','interpreter','latex')
legend('$i = 1$','$i = 2$','$i = 3$','$i = 4$','$i = 5$','$i = 6$',...
    'interpreter','latex','box','off','orientation','horizontal','location','best')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
box on
grid on

plot_ampl_spectrum(Ymag,freq,yL,y0)
