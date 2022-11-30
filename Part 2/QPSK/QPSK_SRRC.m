close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates the BER and the SER numerically for srrc shaped QPSK signals with different SNR values 
% in AWGN channel.
% The program compares the created values with the theoretical results by a graph of
% the BER and the SER as a function of the SNR values in dB.

nos = 1e6; %% number of signals
bss1 = rand(1,nos)>0.5; %% creation of the inphase logical bit stream source randomly
bss2 = rand(1,nos)>0.5; %% creation of the quadrature logical bit stream source randomly


%%creation of the QPSK with seperation of the inphase and the quadrature
I = (2*bss1)-1; 
Q = (2*bss2)-1;
s = I+1j*Q;

%%defining the srrc parameters
span = 20;
sps = 6;
beta = 0.1;
offset = 0;

%%upsampling the given signal and convoluting it with a srrc shaped filter
[pulse_shaped_signal,srrc,tail_elements] = upsample_srrc_conv(s,nos,sps,span,beta);
signal_length = length(pulse_shaped_signal);

SNR_limit_value = 10; %% defining the SNR limit value

gamma_dB = 1:1:SNR_limit_value; %% defining the SNR values in dB
gamma_b = 10.^((gamma_dB)/10);  %% calculation the SNR in linear values
N0 = 1./gamma_b;  %% calculation of the N0 for all SNR values
default_noise = randn(1,nos*sps)+1j*(randn(1,nos*sps)); %% creation of the default noise 

for i = 1:SNR_limit_value
    
    n = sqrt(N0(i)/2)*default_noise; %% calculation of the noise with different N0 values
    r = pulse_shaped_signal + n; %% adding the noise to the symbol
    
    %%convolution of the given signal with a srrc shaped filter and downsampling it
    down_sampled = srrc_conv_downsample(r,srrc,tail_elements,signal_length,sps,offset);

    %%decision for s according to decision law for QPSK
    %%with separate decision for I and Q
    bss1_rx = real(down_sampled)>0;
    bss2_rx = imag(down_sampled)>0;
    
    %%creation of the bit error vectors
    e1 = (bss1~=bss1_rx);
    e2 = (bss2~=bss2_rx);
    
    Nb_error = sum(e1) + sum(e2); %% summing all bit errors in the vector
    BER(i) = Nb_error/(2*nos); %% calculation of the BER 
    
    e_symbol = or(e1,e2); %%  creation of the symbol error vectors
    Ns_error = sum(e_symbol); %% summing all symbol errors in the vector
    SER(i) = Ns_error/nos; %% calculation of the SER
    
end

%%calculation of the theoretical BER and SER
theoretical_BER = qfunc(sqrt(2*gamma_b));
theoretical_SER = 2*qfunc(sqrt(2*gamma_b)) - (qfunc((sqrt(2*gamma_b)))).^2;

%%plotting the graph of the theoretical and the numerical BER and SER as a function of the SNR values in dB.
mod_title = 'QPSK Pulse Shaped Signal';
y_limits = [1e-6 1];
figure(1)
BER_SER_plot(mod_title,gamma_dB,BER,theoretical_BER,SER,theoretical_SER,SNR_limit_value,y_limits);

%%plotting the constellation diagrams
figure(2)
constellation_title = 'QPSK: Pulse Shaped Signal Without Noise';  
display_constellation_QPSK(pulse_shaped_signal,constellation_title)
figure(3)
constellation_title = 'QPSK: Pulse Shaped Signal With Noise';  
display_constellation_QPSK(r,constellation_title)
figure(4)
constellation_title = 'QPSK: Downsampled Signal';  
display_constellation_QPSK(down_sampled,constellation_title)