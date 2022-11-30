clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates the BER numerically for BPSK signals with different SNR values 
% in AWGN channel.
% The program compares the created values with the theoretical results by a graph of
% the BER as a function of the SNR values in dB.

nos = 1e6;  %% number of signals
bss = rand(1,nos)>0.5;  %% creation of the logical bit stream source randomly
s = (2*bss)-1;  %% creation of the BPSK symbols

SNR_limit_value = 10; %% defining the SNR limit value

gamma_dB = 1:1:SNR_limit_value;  %% defining the SNR values in dB
gamma_b = 10.^((gamma_dB)/10);  %% calculation the SNR in linear values
N0 = 1./gamma_b;    %% calculation of the N0 for all SNR values
default_noise = randn(1,nos); %% creation of the default noise 

for i = 1:SNR_limit_value 
    
    n = sqrt(N0(i)/2)*default_noise; %% calculation of the noise with different N0 values
    r = s + n; %% adding the noise to the symbol
    
    bss_rx = r>0; %% decision for s according to decision law for BPSK
    e = (bss~=bss_rx); %% creation of the bit error vector 
    Nb_error = sum(e); %% summing all errors in the vector
    BER(i) = Nb_error/nos; %% calculation of the BER
    
end

theoretical_BER = qfunc(sqrt(2*gamma_b));%% calculation of the theoretical BER 

%%plotting the graph of the theoretical and the numerical BER as a function of the SNR values in dB.
figure(1)
BER_plot = semilogy(gamma_dB,BER,'g',gamma_dB,theoretical_BER,'--ob');
BER_plot(1).LineWidth = 2;
BER_plot(2).LineWidth = 2;
xlim([1 10])
ylim([1e-6 1e-1])
xlabel('SNR [dB]')
ylabel('BER')
title('BPSK in awgn channel BER')
legend('Numeric BER','Theoretical BER')
grid on
