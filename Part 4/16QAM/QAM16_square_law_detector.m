close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates srrc shaped 16QAM signals and transmits them through 
% an optical fiber with dispersion and chromatic dispersion compensation in the DSP level.
% The program then calculates the BER and the SER numerically with different SNR values in AWGN channel.
% The program compares the created values with the theoretical results by a graph of
% the BER and the SER as a function of the SNR values in dB.

nos = 1e6; %% number of signals

%%creation of the 16QAM symbols using a designated function
[s,bssI_1,bssI_2,bssQ_1,bssQ_2] = define16qam(nos);

%%defining the srrc parameters
span = 20;
sps = 6;
beta = 0.1;
offset = 0;

%%upsampling the given signal and convoluting it with a srrc shaped filter
[pulse_shaped_signal,srrc,tail_elements] = upsample_srrc_conv(s,nos,sps,span,beta);

%%definition of the symbol rate and sampling frequencies for the optical fiber
signal_length = length(pulse_shaped_signal);
SR = 30e9; % [baud]
Fs = SR*sps;
Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));

%%definition of the sampling frequencies for the CDC
Fs2 = SR*2;
Fv2 = Fs2*(-0.5:(1/(2*nos)):(0.5-1/(2*nos)));

%%definition of the fiber parameters
D1 = 17e-15; %% [s/nm-m]
L1 = 0; %% [m]
D2 = -100e-15; %%[s/nm-m]
L2 = 0; %% [m]

%%transmission of the signal through the optical fiber
y = fiber(pulse_shaped_signal,Fv,D1,L1);

%%convertion of the signal to SSB configuration with carrier when CSPR=10
[ssbshape,carrier] = ssbshaper(y,beta,SR,Fs);

%%plotting the fourier transform of each part on a logarithmic scale
figure(1)
title1 = 'Es (f)';
fft_plot_log(y,Fs,title1);
figure(2)
title2 = 'Eo (f + B/2)';
fft_plot_log(carrier,Fs,title2);
figure(3)
title3 = 'Es (f)  +  Eo (f + B/2)';
fft_plot_log(ssbshape,Fs,title3);

M = 16;      

SNR_limit_value = 10; %% defining the SNR limit value
gamma_dB = SNR_limit_value; %% defining the SNR values in dB

%%calculation of the N0 for all SNR values and creation of the default noise 
[N0,default_noise,gamma_b] = create_N0_defaultnoise(s,M,nos,sps,gamma_dB);

n = sqrt(N0/2)*default_noise; %% calculation of the noise with different N0 values

r = ssbshape + n; %% adding the noise to the symbol

%%using squarelaw detection to detect the signal and creating each part of the result separately
squarelaw = (abs(ssbshape)).^2;
squarelaw_signal = 2*real(y.*conj(carrier));
squarelaw_carrier = (abs(carrier)).^2;
squarelaw_SSBI = (abs(y)).^2;

%%plotting the fourier transform of each part on a logarithmic scale
%%we present the graphs in frequency domain without the noise
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

%%convolution of the given signal with a srrc shaped filter and downsampling it
down_sampled = srrc_conv_downsample(squarelaw_with_noise,srrc,tail_elements,signal_length,3,offset);

%%chromatic dispersion compensation in the dsp level
y_CDC = fiber(down_sampled,Fv2,D2,L2);

%%downsampling to the digital level
down_sampled_2 = downsample(y_CDC,2);

%%estimating the inphase and quadrature bits from the noisy signal using a designated function
[bssI_1_rx,bssI_2_rx,bssQ_1_rx,bssQ_2_rx] = decision16qam(down_sampled_2);

%%creation of the bit error vectors
eI1 = (bssI_1~=bssI_1_rx);
eI2 = (bssI_2~=bssI_2_rx);

eQ1 = (bssQ_1~=bssQ_1_rx);
eQ2 = (bssQ_2~=bssQ_2_rx);

Nb_error = sum(eI1) + sum(eI2) + sum(eQ1) + sum(eQ2); %% summing all bit errors in the vector
BER = Nb_error/(4*nos); %% calculation of the BER 

e_symbol = eI1 | eI2 | eQ1 | eQ2; %%  creation of the symbol error vectors
Ns_error = sum(e_symbol); %% summing all symbol errors in the vector
SER = Ns_error/nos; %% calculation of the SER

%%calculation of the theoretical BER and SER
theoretical_SER =(1-(1-2*(sqrt(M)-1)/sqrt(M)*qfunc(sqrt(3*log2(M)/(M-1)*gamma_b))).^2);
theoretical_BER = theoretical_SER/log2(M);
