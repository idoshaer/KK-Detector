close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program creates a srrc shaped 64-QAM signal and transmits it through an AWGN channel.  
% The signal is detected by using the squarelaw.
% The program implements the KK algorithm to cancel the SSBI for
% recreating the signal's spectrum.
% The program then calculates the BER numerically and shows the BER as the function
% of the total SNR for different CSPR values.


nos = 1e6; %% number of signals

M = 64; %number of symbols (64-QAM constellation)      

%%creating the 64QAM bit streams and defining the symbols using a
%%designated function
[s,bssI_1,bssI_2,bssI_3,bssQ_1,bssQ_2,bssQ_3] = define64qam(nos);

%%defining the srrc parameters
span = 20;
initial_sps = 1;
beta = 0.1;

%%definition of the fiber parameters
D1 = 17e-15; %% [s/nm-m]
L1 = 100000; %% [m]
D2 = -100e-15; %%[s/nm-m]
L2 = 17e3; %% [m]

%%upsampling the signal to sps=6
[upsampled_signal,sps] = upsample_sps(s,initial_sps,6);

%%convoluting the signal with a srrc shaped filter
[pulse_shaped_signal,srrc1] = srrc_conv(upsampled_signal,sps,span,beta);

%%definition of the symbol rate and sampling frequencies for the optical fiber
signal_length = length(pulse_shaped_signal);
SR = 30e9; % [baud]
Fs = SR*sps;
Fv = Fs*(-0.5:(1/signal_length):(0.5-1/signal_length));

%%transmission of the signal through the first optical fiber
signal_fiber = fiber(pulse_shaped_signal,Fv,D1,L1);

%%definition of the sampling frequencies for the CDC
Fs2 = SR*2;
Fv2 = Fs2*(-0.5:(1/(2*nos)):(0.5-1/(2*nos)));

gamma_dB = 1:2:25;

gamma_dB_total_matrix = [];

%%defining the sps values before KK receiver
desired_sps = 3.5;

%%defining the CSPR values
CSPR_dB = [5 8 10 13 15];%% in dB

BER_matrix = [];

for k = gamma_dB
    
    %%calculation of the N0 for all SNR values and creation of the default noise 
    [N0,default_noise,gamma_b] = create_N0_defaultnoise(s,M,nos,sps,k);

    %%calculation of the noise with different N0 values
    n = sqrt(N0/2)*default_noise; 
    
    %%defining the SNR in dB
    SNR_signal_lin = (db2mag(k))^2;
    
    gamma_dB_total_vector = [];
    BER_vector = [];

    for i = CSPR_dB

        CSPR_lin = (db2mag(i))^2;
        gamma_dB_total = 0.5*mag2db(SNR_signal_lin*(1+CSPR_lin)); 
        gamma_dB_total_vector = [gamma_dB_total_vector,gamma_dB_total];

        %%converting the signal to SSB configuration with carrier for different values of CSPR
        [ssbshape,carrier] = ssbshaper(signal_fiber,beta,SR,Fs,i);

        %%adding the noise to the symbol
        r = ssbshape + n; 

        %%squarelaw detection
        squarelaw_with_noise = (abs(r)).^2;

        %%resampling the signal to different sps values
        [P,Q] = numden(sym(desired_sps/6));
        resample1 = resample(squarelaw_with_noise,double(P),double(Q),2);
        sps = sps*double(P)/double(Q);
        Fs_1 = SR*sps;

        %%implementing the KK algorithm
        shifted_frequency = KK_algorithm(resample1,Fs_1,beta,SR);

        %%resampling the signal to sps=2 
        [P2,Q2] = numden(sym(1/desired_sps));
        resample2 = resample(shifted_frequency,2*double(P2),double(Q2),2);
        sps = 2;

        %%chromatic dispersion compensation in the dsp level
        signal_CDC = fiber(resample2,Fv2,D2,L2); 

        %%convoluting with a srrc shaped filter
        [pulse_shaped_signal_new,srrc2] = srrc_conv(signal_CDC,2,span,beta);

        %%normalization of the signal after matched filter
        conv_srrc = conv(srrc1,srrc2);
        pulse_shaped_signal_new = pulse_shaped_signal_new/max(conv_srrc);

        %%downsampling to the digital level
        [downsampled_signal_3,sps] = downsample_sps(pulse_shaped_signal_new,sps,2,0);

        %%estimating the inphase and quadrature bits using a designated function
        [bssI_1_rx,bssI_2_rx,bssI_3_rx,bssQ_1_rx,bssQ_2_rx,bssQ_3_rx] = decision64qam(downsampled_signal_3);

        %%creation of the bit error vectors
        eI1 = (bssI_1~=bssI_1_rx);
        eI2 = (bssI_2~=bssI_2_rx);
        eI3 = (bssI_3~=bssI_3_rx);

        eQ1 = (bssQ_1~=bssQ_1_rx);
        eQ2 = (bssQ_2~=bssQ_2_rx);
        eQ3 = (bssQ_3~=bssQ_3_rx);

        Nb_error = sum(eI1) + sum(eI2) + sum(eI3) + sum(eQ1) + sum(eQ2)+ sum(eQ3); %% summing all bit errors in the vector
        BER = Nb_error/(6*nos); %% calculation of the BER 
        BER_vector = [BER_vector,BER];

        %%resetting the value of the sps to 6 
        sps = 6;

    end
    
    gamma_dB_total_matrix = [gamma_dB_total_matrix;gamma_dB_total_vector];
    BER_matrix = [BER_matrix;BER_vector];

end    


figure(1)
BER_plot = semilogy(gamma_dB_total_matrix(:,1),BER_matrix(:,1),'-b');
xlim([5 45]);
ylim([4e-3 4e-1]);
xlabel('SNR total [dB]');
ylabel('BER');
title('64QAM: BER as a function of the total SNR for different CSPR values when sps dsp = 3.5')
grid on;

hold on
BER_plot2 = semilogy(gamma_dB_total_matrix(:,2),BER_matrix(:,2),'-r');
BER_plot3 = semilogy(gamma_dB_total_matrix(:,3),BER_matrix(:,3),'-y');
BER_plot4 = semilogy(gamma_dB_total_matrix(:,4),BER_matrix(:,4),'-k');
BER_plot5 = semilogy(gamma_dB_total_matrix(:,5),BER_matrix(:,5),'-g');

legend('CSPR=5dB','CSPR=8dB','CSPR=10dB','CSPR=13dB','CSPR=15dB');
