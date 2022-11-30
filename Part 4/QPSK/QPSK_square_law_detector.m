close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates srrc shaped QPSK signals and converts the signal to SSB configuration.
% The program then uses squarelaw detection to detect the signal.
% The program presents each step in frequency domain.
% Finally, the program calculates the BER and the SER numerically in AWGN channel.

nos = 1e6; %% number of signals
bss1 = rand(1,nos)>0.5; %% creation of the inphase logical bit stream source randomly
bss2 = rand(1,nos)>0.5; %% creation of the quadrature logical bit stream source randomly


%%% creation of the QPSK with seperation of the inphase and the quadrature
I = (2*bss1)-1; 
Q = (2*bss2)-1;
s = I+1j*Q;

%%% defining the srrc parameters
span = 20;
sps = 6;
beta = 0.1;
offset = 0;

%%% upsampling the given signal and convoluting it with a srrc shaped filter
[pulse_shaped_signal,srrc,tail_elements] = upsample_srrc_conv(s,nos,sps,span,beta);

%%% definition of the different parameters
signal_length = length(pulse_shaped_signal);
SR = 30e9; % [baud]
Fs = SR*sps;
Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));

%%% definition of the sampling frequencies for the CDC
Fs2 = SR*2;
Fv2 = Fs2*(-0.5:(1/(2*nos)):(0.5-1/(2*nos)));

%%% definition of the fiber parameters
D1 = 17e-15; %% [s/nm-m]
L1 = 0; %% [m]
D2 = -100e-15; %%[s/nm-m]
L2 = 0; %% [m]

%%% transmission of the signal through the optical fiber
y = fiber(pulse_shaped_signal,Fv,D1,L1);

%%% convertion of the signal to SSB configuration with carrier when CSPR=10
[ssbshape,carrier] = ssbshaper(y,beta,SR,Fs);

%%% plotting the fourier transform of each part on a logarithmic scale
figure(1)
title1 = 'Es (f)';
fft_plot_log(y,Fs,title1);
figure(2)
title2 = 'Eo (f + B/2)';
fft_plot_log(carrier,Fs,title2);
figure(3)
title3 = 'Es (f)  +  Eo (f + B/2)';
fft_plot_log(ssbshape,Fs,title3);

SNR = 10; %% defining the SNR limit value
gamma_b = 10^((SNR)/10);  %% calculation the SNR in linear values
N0 = 1/gamma_b;  %% calculation of the N0 for all SNR values
default_noise = randn(1,signal_length)+1j*(randn(1,signal_length)); %% creation of the default noise 
n = sqrt(N0/2)*default_noise; %% calculation of the noise with different N0 values
r = ssbshape + n; %% adding the noise to the symbol

%%% using squarelaw detection to detect the signal and creating each part of the result separately
squarelaw = (abs(ssbshape)).^2;
squarelaw_signal = 2*real(y.*conj(carrier));
squarelaw_carrier = (abs(carrier)).^2;
squarelaw_SSBI = (abs(y)).^2;

%%% plotting the fourier transform of each part on a logarithmic scale
%%% we present the graphs in frequency domain without the noise
figure(4)
title4 = 'Fourier Transform of  | E (t) | ^2';
fft_plot_log(squarelaw,Fs,title4);
figure(41)
title41 = 'Fourier Transform of  2 * Re [Es (t) * Eo * exp(j*pi*B*t)]';
fft_plot_log(squarelaw_signal,Fs,title41);
figure(42)
title42 = 'Fourier Transform of  | Eo (t) | ^2';
fft_plot_log(squarelaw_carrier,Fs,title42);
figure(43)
title43 = 'Fourier Transform of  | Es (t) | ^2';
fft_plot_log(squarelaw_SSBI,Fs,title43);

squarelaw_with_noise = (abs(r)).^2;

%%% convolution of the given signal with a srrc shaped filter and downsampling it
down_sampled = srrc_conv_downsample(squarelaw_with_noise,srrc,tail_elements,signal_length,3,offset);

%%% chromatic dispersion compensation in the dsp level
y_CDC = fiber(down_sampled,Fv2,D2,L2);

%%% downsampling to the digital level
down_sampled_2 = downsample(y_CDC,2);

%%% decision for s according to decision law for QPSK
%%% with separate decision for I and Q
bss1_rx = real(down_sampled_2)>0;
bss2_rx = imag(down_sampled_2)>0;

%%% creation of the bit error vectors
e1 = (bss1~=bss1_rx);
e2 = (bss2~=bss2_rx);

Nb_error = sum(e1) + sum(e2); %% summing all bit errors in the vector
BER = Nb_error/(2*nos); %% calculation of the BER 

e_symbol = or(e1,e2); %%  creation of the symbol error vectors
Ns_error = sum(e_symbol); %% summing all symbol errors in the vector
SER = Ns_error/nos; %% calculation of the SER

%%% calculation of the theoretical BER and SER
theoretical_BER = qfunc(sqrt(2*gamma_b));
theoretical_SER = 2*qfunc(sqrt(2*gamma_b)) - (qfunc((sqrt(2*gamma_b)))).^2;
