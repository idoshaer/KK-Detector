function [] = fft_plot_log(time_domain_signal,Fs,graph_title)
%%this function performs fourier transform to a signal in time domain and
%%plots it on a logarithmic scale

signal_length = length(time_domain_signal);
Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));
signal_fft = fftshift(abs(fft(time_domain_signal)/signal_length));
plot(Fv,10*log10(signal_fft));
xlabel('frequency [Hz]')
ylabel(graph_title)
title(graph_title)
grid on


end

