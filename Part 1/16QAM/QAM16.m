close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates the BER and the SER numerically for 16QAM signals with different SNR values 
% in AWGN channel.
% The program compares the created values with the theoretical results by a graph of
% the BER and the SER as a function of the SNR values in dB.

nos = 1e6; %% number of signals

%%creation of the 16QAM symbols using a designated function
[s,bssI_1,bssI_2,bssQ_1,bssQ_2] = define16qam(nos);

M = 16;

SNR_limit_value = 10; %% defining the SNR limit value

%%calculation of the N0 for all SNR values and creation of the default noise 
[N0,default_noise,gamma_b,gamma_dB] = create_N0_defaultnoise(s,M,nos,1,SNR_limit_value);


for i = 1:SNR_limit_value
    
    n = sqrt(N0(i)/2)*default_noise; %% calculation of the noise with different N0 values
    r = s + n; %% adding the noise to the symbol
    
    %%estimating the inphase and quadrature bits from the noisy signal using a designated function
    [bssI_1_rx,bssI_2_rx,bssQ_1_rx,bssQ_2_rx] = decision16qam(r);

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
mod_title = '16QAM in AWGN Channel';
y_limits = [1e-4 1];
figure(1)
BER_SER_plot(mod_title,gamma_dB,BER,theoretical_BER,SER,theoretical_SER,SNR_limit_value,y_limits);

%%plotting the constellation diagrams with and without noise
figure(2)
constellation_title = '16QAM: Without Noise';  
display_constellation_16QAM(s,constellation_title)
figure(3)
constellation_title = '16QAM: With Noise';  
display_constellation_16QAM(r,constellation_title)
