close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates srrc shaped 16QAM signals and transmits them through an optical fiber with dispersion and DCF. 
% The program then calculates the BER for different length deviations of
% the optical fiber in 3 different symbol rates.

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

signal_length = length(pulse_shaped_signal);

%%definition of the fiber parameters
D1 = 17e-15; %% [s/nm-m]
L1 = 100000; %% [m]
delta_L = (-1e3:25:1e3); %% [m]

%%calculation of the average energy per symbol and bit
M = 16;

SNR_limit_value = 10; %% defining the SNR limit value
gamma_dB = SNR_limit_value; %% defining the SNR values in dB

%%calculation of the N0 for all SNR values and creation of the default noise 
[N0,default_noise,gamma_b] = create_N0_defaultnoise(s,M,nos,sps,gamma_dB);

n = sqrt(N0/2)*default_noise; %% calculation of the noise

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

        %%estimating the inphase and quadrature bits from the noisy signal using a designated function
        [bssI_1_rx,bssI_2_rx,bssQ_1_rx,bssQ_2_rx] = decision16qam(down_sampled);
    
        %%creation of the bit error vectors
        eI1 = (bssI_1~=bssI_1_rx);
        eI2 = (bssI_2~=bssI_2_rx);

        eQ1 = (bssQ_1~=bssQ_1_rx);
        eQ2 = (bssQ_2~=bssQ_2_rx);
    
        Nb_error = sum(eI1) + sum(eI2) + sum(eQ1) + sum(eQ2); %% summing all bit errors in the vector
        BER = Nb_error/(4*nos); %% calculation of the BER 
        BER_vector = [BER_vector,BER]; %% allocating the BER values to different length deviations
        
    end
    
    BER_matrix = [BER_matrix;BER_vector]; %% creation of the BER vectors for different symbol rates
    
end

mod_title = '16QAM';

%%plotting the BER for different length deviations of the optical fiber in 3 different symbol rates
BER_delta_L_SR_plot(delta_L,BER_matrix,mod_title);
