function [magnitude,phase,f,x] = FFT(x0,nL,L,Fs,opts)

% Function to compute the FFT of a signal x0
% Written: Dinesh Krishnamoorthy , Aug 2020; dineshk@ntnu.no

if nargin<5
    opts.plot = 0;
end

x = x0 - movingavg(x0,L/2,0); % detrend using moving average
NFFT = 2^nextpow2(nL*L); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);

X = fft(x,NFFT)/NFFT;

magnitude = 2*abs(X(1:NFFT/2+1));
phase = angle(X(1:NFFT/2+1))*180/pi;

if opts.plot
    % Plot single-sided amplitude spectrum
    figure(12)
    subplot(212)
    stem(f,magnitude,'linewidth',2)
    title('Single-Sided Amplitude Spectrum of $J_0(t)$','Interpreter','Latex')
    xlabel('Frequency (Hz)','Interpreter','Latex')
    ylabel('$2|\mathcal{J}(\omega)|$','Interpreter','Latex')
    xlim([0,200])
    grid on
    axs = gca;
    axs.FontSize = 14;
    % axs.FontName = 'SanSerif';
    axs.TickLabelInterpreter = 'latex';
    yticks([0,0.2,0.4,0.6,0.8,1])
    
    subplot(211)
    plot(x,'k','linewidth',1.75)
    xlabel('Time $t$ ','Interpreter','Latex')
    ylabel('$J_0$','Interpreter','Latex')
    axs = gca;
    axs.FontSize = 14;
    % axs.FontName = 'SanSerif';
    axs.TickLabelInterpreter = 'latex';
    
end
end
