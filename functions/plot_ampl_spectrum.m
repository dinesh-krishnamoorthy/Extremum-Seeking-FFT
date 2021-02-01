function plot_ampl_spectrum(Xmag,f,x,x0)

% Function to plot the single sided amplitude spectrum
% Written: Dinesh Krishnamoorthy Aug 2020

figure()
clf
subplot(311)
plot(x,'k','linewidth',1.75)
xlabel('Time $t$ ','Interpreter','Latex')
ylabel('$J$','Interpreter','Latex')
axs = gca;
axs.FontSize = 14;
% axs.FontName = 'SanSerif';
axs.TickLabelInterpreter = 'latex';
grid on

subplot(312)
plot(x0,'k','linewidth',1.75)
xlabel('Time $t$ ','Interpreter','Latex')
ylabel('$J_0$','Interpreter','Latex')
axs = gca;
axs.FontSize = 14;
% axs.FontName = 'SanSerif';
axs.TickLabelInterpreter = 'latex';
grid on

subplot(313)
stem(f,Xmag,'linewidth',2) 
title('Single-Sided Amplitude Spectrum of $J_0(t)$','Interpreter','Latex')
xlabel('Frequency (Hz)','Interpreter','Latex')
ylabel('$2|\mathcal{J}(\omega)|$','Interpreter','Latex')
grid on
axs = gca;
axs.FontSize = 14;
% axs.FontName = 'SanSerif';
axs.TickLabelInterpreter = 'latex';
% yticks([0,0.2,0.4,0.6,0.8,1])

