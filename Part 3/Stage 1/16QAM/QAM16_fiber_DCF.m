close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates srrc shaped 16QAM signals and transmits them through an optical fiber with dispersion and DCF.
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

D1 = 17e-15; %% [s/nm-m]
L1 = 100000; %% [m]

%%transmission of the signal through the first optical fiber
y = fiber(pulse_shaped_signal,Fv,D1,L1);

%%definition of the DCF parameters
D2 = -100e-15; %%[s/nm-m]
L2 = 17e3; %% [m]

%%transmission of the signal through the DCF
y_DCF = fiber(y,Fv,D2,L2);

M = 16;      

SNR_limit_value = 10; %% defining the SNR limit value
gamma_dB = 1:1:SNR_limit_value; %% defining the SNR values in dB

%%calculation of the N0 for all SNR values and creation of the default noise 
[N0,default_noise,gamma_b] = create_N0_defaultnoise(s,M,nos,sps,gamma_dB);

for i = 1:SNR_limit_value
    
    n = sqrt(N0(i)/2)*default_noise; %% calculation of the noise with different N0 values
    r = y_DCF + n; %% adding the noise to the symbol
    
    %%convolution of the given signal with a srrc shaped filter and
    %%downsampling it
    down_sampled = srrc_conv_downsample(r,srrc,tail_elements,signal_length,sps,offset);
    
    %%estimating the inphase and quadrature bits from the noisy signal using a designated function
    [bssI_1_rx,bssI_2_rx,bssQ_1_rx,bssQ_2_rx] = decision16qam(down_sampled);
    
    %%creation of the bit error vectors
    eI1 = (bssI_1~=bssI_1_rx);
    eI2 = (bssI_2~=bssI_2_rx);
    
    eQ1 = (bssQ_1~=bssQ_1_rx);
    eQ2 = (bssQ_2~=bssQ_2_rx);
    
 
    Nb_error = sum(eI1) + sum(eI2) + sum(eQ1) + sum(eQ2); %% summing all bit errors in the vector
    BER(i) = Nb_error/(4*nos); %% calculation of the BER 
    
    e_symbol = eI1 | eI2 | eQ1 | eQ2; %%  creation of the symbol error vectors
    Ns_error = sum(e_symbol); %% summing all symbol errors in the vector
    SER(i) = Ns_error/nos; %% calculation of the SER

end

%%calculation of the theoretical BER and SER
theoretical_SER =(1-(1-2*(sqrt(M)-1)/sqrt(M)*qfunc(sqrt(3*log2(M)/(M-1)*gamma_b))).^2);
theoretical_BER = theoretical_SER/log2(M);

%%plotting the graph of the theoretical and the numerical BER and SER as a function of the SNR values in dB.
mod_title = '16QAM Through Fiber and DCF';
y_limits = [1e-4 1];
figure(1)
BER_SER_plot(mod_title,gamma_dB,BER,theoretical_BER,SER,theoretical_SER,SNR_limit_value,y_limits);

%%plotting the constellation diagrams
figure(2)
constellation_title = '16QAM: Pulse Shaped Signal Without Noise';  
display_constellation_16QAM(pulse_shaped_signal,constellation_title)
figure(3)
constellation_title = '16QAM: Transmission of the Signal Through the First Optical Fiber';  
display_constellation_16QAM(y,constellation_title)
figure(4)
constellation_title = '16QAM: Transmission of the Signal Through the DCF';  
display_constellation_16QAM(y_DCF,constellation_title)
