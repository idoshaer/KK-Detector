close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates a srrc shaped 16-QAM signal and transmits it through an AWGN channel.  
% The signal is detected by using the squarelaw.
% The program implements the KK algorithm to cancel the SSBI for
% recreating the signal's spectrum.
% The program then calculates the BER numerically and shows the spectrum graphs
% of different stages in the process and the constellation diagram before
% the decision.


nos = 1e6; %% number of signals

%%creation of the 16QAM symbols using a designated function
[s,bssI_1,bssI_2,bssQ_1,bssQ_2] = define16qam(nos);

%%defining the srrc parameters
span = 20;
initial_sps = 1;
beta = 0.1;

%%upsampling the given signal
interpolation_index = 6;
[upsampled_signal,sps] = upsample_sps(s,initial_sps,interpolation_index);

%%convoluting the signal with a srrc shaped filter
[pulse_shaped_signal,srrc1] = srrc_conv(upsampled_signal,sps,span,beta);

%%definition of the symbol rate and sampling frequencies for the optical fiber
signal_length = length(pulse_shaped_signal);
SR = 30e9; % [baud]
Fs = SR*sps;
Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));

%%convertion of the signal to SSB configuration with carrier when CSPR=10
CSPR_dB = 10; %% in dB
[ssbshape,carrier] = ssbshaper(pulse_shaped_signal,beta,SR,Fs,CSPR_dB);

%%plotting the fourier transform on a logarithmic scale
figure(1)
title1 = 'Es (f)  +  Eo (f + B/2)';
fft_plot_log(ssbshape,Fs,title1);

M = 16; %number of symbols (16-QAM constellation)      

SNR_limit_value = 10; %% defining the SNR limit value
gamma_dB = SNR_limit_value; %% defining the SNR values in dB

%%calculation of the N0 for all SNR values and creation of the default noise 
[N0,default_noise,gamma_b] = create_N0_defaultnoise(s,M,nos,sps,gamma_dB);

n = sqrt(N0/2)*default_noise; %% calculation of the noise with different N0 values

r = ssbshape + n; %% adding the noise to the symbol

%%squarelaw detection
squarelaw_with_noise = (abs(r)).^2;
figure(2)
title2 = 'Fourier Transform of  | E (t) | ^2';
fft_plot_log(squarelaw_with_noise,Fs,title2);

%%downsampling to sps = 3
decimation_index = 2;
offset = 0;
[downsampled_signal,sps] = downsample_sps(squarelaw_with_noise,sps,decimation_index,offset);
figure(3)
title3 = 'downsampled 1';
fft_plot_log(downsampled_signal,Fs,title3);

%%implementing the KK algorithm
shifted_frequency = KK_algorithm(downsampled_signal,Fs,beta,SR);

%%convoluting with a srrc shaped filter
[pulse_shaped_signal_new,srrc2] = srrc_conv(shifted_frequency,sps,span,beta);

%%normalization of the signal after matched filter
conv_srrc = conv(srrc1,srrc2);
pulse_shaped_signal_new = pulse_shaped_signal_new/max(conv_srrc);

%%downsampling to the digital level
decimation_index = 3;
[downsampled_signal_2,sps] = downsample_sps(pulse_shaped_signal_new,sps,decimation_index,offset);
figure(6)
title6 = 'downsampled signal 2';
fft_plot_log(downsampled_signal_2,Fs,title6);

%%estimating the inphase and quadrature bits from the noisy signal using a designated function
[bssI_1_rx,bssI_2_rx,bssQ_1_rx,bssQ_2_rx] = decision16qam(downsampled_signal_2);

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

%%plotting the constellation diagram
figure(7)
constellation_title = '16QAM: After KK Algorithm (CSPR=10dB SNR=10dB)';  
display_constellation_16QAM(downsampled_signal_2,constellation_title)

