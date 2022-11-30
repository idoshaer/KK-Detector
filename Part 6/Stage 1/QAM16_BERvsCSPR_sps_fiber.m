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

%%definition of the fiber parameters
D1 = 17e-15; %% [s/nm-m]
L1 = 100000; %% [m]
D2 = -100e-15; %%[s/nm-m]
L2 = 17e3; %% [m]

M = 16; %number of symbols (16-QAM constellation)      

gamma_dB = 15; %% defining the SNR values in dB

%%upsampling the given signal
interpolation_index = [4 6 8 10 12];

BER_matrix = [];


for k = interpolation_index
    
    [upsampled_signal,sps] = upsample_sps(s,initial_sps,k);

    %%convoluting the signal with a srrc shaped filter
    [pulse_shaped_signal,srrc1] = srrc_conv(upsampled_signal,sps,span,beta);
    
    %%definition of the symbol rate and sampling frequencies for the optical fiber
    signal_length = length(pulse_shaped_signal);
    SR = 30e9; % [baud]
    Fs = SR*sps;
    Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));
    
    %%transmission of the signal through the first optical fiber
    signal_fiber = fiber(pulse_shaped_signal,Fv,D1,L1);
    
    %%definition of the decimation parameters
    decimation_index_1 = 2;
    decimation_index_2 = k/decimation_index_1;
    offset = 0;
    
    %%definition of the sampling frequencies for the CDC
    Fs2 = SR*decimation_index_2;
    Fv2 = Fs2*(-0.5:(1/(decimation_index_2*nos)):(0.5-1/(decimation_index_2*nos)));

    %%calculation of the N0 for all SNR values and creation of the default noise 
    [N0,default_noise,gamma_b] = create_N0_defaultnoise(s,M,nos,sps,gamma_dB);

    n = sqrt(N0/2)*default_noise; %% calculation of the noise with different N0 values

    %%convertion of the signal to SSB configuration with carrier when CSPR=10
    CSPR_dB = 1:2:25;%% in dB

    BER_vector = [];

    for i = CSPR_dB

        [ssbshape,carrier] = ssbshaper(signal_fiber,beta,SR,Fs,i);

        r = ssbshape + n; %% adding the noise to the symbol

        %%squarelaw detection
        squarelaw_with_noise = (abs(r)).^2;

        %%downsampling according to the first decimation index
        [downsampled_signal,sps] = downsample_sps(squarelaw_with_noise,sps,decimation_index_1,offset);

        %%implementing the KK algorithm
        shifted_frequency = KK_algorithm(downsampled_signal,Fs/2,beta,SR);
        
        %%chromatic dispersion compensation in the dsp level
        signal_CDC = fiber(shifted_frequency,Fv2,D2,L2); 
        
        %%convoluting with a srrc shaped filter
        [pulse_shaped_signal_new,srrc2] = srrc_conv(signal_CDC,sps,span,beta);

        %%normalization of the signal after matched filter
        conv_srrc = conv(srrc1,srrc2);
        pulse_shaped_signal_new = pulse_shaped_signal_new/max(conv_srrc);

        %%downsampling to the digital level
        [downsampled_signal_2,sps] = downsample_sps(pulse_shaped_signal_new,sps,decimation_index_2,offset);

        %%estimating the inphase and quadrature bits from the noisy signal using a designated function
        [bssI_1_rx,bssI_2_rx,bssQ_1_rx,bssQ_2_rx] = decision16qam(downsampled_signal_2);

        %%creation of the bit error vectors
        eI1 = (bssI_1~=bssI_1_rx);
        eI2 = (bssI_2~=bssI_2_rx);

        eQ1 = (bssQ_1~=bssQ_1_rx);
        eQ2 = (bssQ_2~=bssQ_2_rx);

        Nb_error = sum(eI1) + sum(eI2) + sum(eQ1) + sum(eQ2); %% summing all bit errors in the vector
        BER = Nb_error/(4*nos); %% calculation of the BER 
        BER_vector = [BER_vector,BER];
        sps = k;

    end
    
    BER_matrix = [BER_matrix;BER_vector];
    
end

figure(1)
BER_plot = semilogy(CSPR_dB,BER_matrix(1,:),'-og');
xlim([0 26]);
ylim([2e-3 4e-1]);
xlabel('CSPR [dB]');
ylabel('BER');
title('BER as a function of CSPR for different sps values - through optical fiber')
grid on;

hold on
BER_plot2 = semilogy(CSPR_dB,BER_matrix(2,:),'-ob');
BER_plot3 = semilogy(CSPR_dB,BER_matrix(3,:),'-oy');
BER_plot4 = semilogy(CSPR_dB,BER_matrix(4,:),'--or');
BER_plot5 = semilogy(CSPR_dB,BER_matrix(5,:),'-ok');

legend('2 dsp sps','3 dsp sps','4 dsp sps','5 dsp sps','6 dsp sps');

