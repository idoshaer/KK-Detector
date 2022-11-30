close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates srrc shaped QPSK signals and transmits them through an optical fiber with dispersion and DCF. 
% The program then calculates the BER for different length deviations of
% the optical fiber in 3 different symbol rates.

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

%%definition of the fiber parameters
D1 = 17e-15; %% [s/nm-m]
L1 = 100000; %% [m]
delta_L = (-5e3:100:5e3); %% [m]

gamma_dB = 7; %% defining the SNR value in dB
gamma_b = 10^((gamma_dB)/10);  %% calculation the SNR in linear value
N0 = 1/gamma_b;  %% calculation of the N0 
n =  sqrt(N0/2)*(randn(1,signal_length)+1j*(randn(1,signal_length))); %% creation of the noise 

BER_matrix = [];

for SR = [15e9 30e9 60e9]
    
    %%creation of the frequency axis
    Fs = SR*sps;
    Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));
    BER_vector = [];
    
    for i = delta_L
        
        %%creation of the length deviations and transmission through the fiber
        L1_v = L1 + i;
        y = fiber(pulse_shaped_signal,Fv,D1,L1_v);
        
        %%definition of the DCF parameters and transmission through the DCF
        D2 = -100e-15; %%[s/nm-m]
        L2 = 17e3; %% [m]
        y_DCF = fiber(y,Fv,D2,L2);
        
        r = y_DCF + n; %% adding the noise to the symbol

        %%convolution of the given signal with a srrc shaped filter and downsampling it
        down_sampled = srrc_conv_downsample(r,srrc,tail_elements,signal_length,sps,offset);

        %%decision for s according to decision law for QPSK with separate decision for I and Q
        bss1_rx = real(down_sampled)>0;
        bss2_rx = imag(down_sampled)>0;

        %%creation of the bit error vectors
        e1 = (bss1~=bss1_rx);
        e2 = (bss2~=bss2_rx);

        Nb_error = sum(e1) + sum(e2); %% summing all bit errors in the vector
        BER = Nb_error/(2*nos); %% calculation of the BER 
        BER_vector = [BER_vector,BER]; %% allocating the BER values to different length deviations
    end
    BER_matrix = [BER_matrix;BER_vector]; %% creation of the BER vectors for different symbol rates
end

%%plotting the BER for different length deviations of the optical fiber in 3 different symbol rates
mod_title = 'QPSK';
BER_delta_L_SR_plot(delta_L,BER_matrix,mod_title);
