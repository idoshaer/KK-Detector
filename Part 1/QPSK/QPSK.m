close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates the BER and the SER numerically for QPSK signals with different SNR values 
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

SNR_limit_value = 10; %% defining the SNR limit value

gamma_dB = 1:1:SNR_limit_value; %% defining the SNR values in dB
gamma_b = 10.^((gamma_dB)/10);  %% calculation the SNR in linear values
N0 = 1./gamma_b;  %% calculation of the N0 for all SNR values
default_noise = randn(1,nos)+1j*(randn(1,nos)); %% creation of the default noise 

for i = 1:SNR_limit_value
    
    n = sqrt(N0(i)/2)*default_noise; %% calculation of the noise with different N0 values
    r = s + n; %% adding the noise to the symbol
    
    %%decision for s according to decision law for QPSK with separate decision for I and Q
    bss1_rx = real(r)>0;
    bss2_rx = imag(r)>0;
    
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

%%plotting the graph of the theoretical and the numerical BER and SER as a
%%function of the SNR values in dB.
mod_title = 'QPSK in AWGN Channel';
y_limits = [1e-6 1];
figure(1)
BER_SER_plot(mod_title,gamma_dB,BER,theoretical_BER,SER,theoretical_SER,SNR_limit_value,y_limits);

%%plotting the constellation diagrams with and without noise
figure(2)
constellation_title = 'QPSK: Without Noise';  
display_constellation_QPSK(s,constellation_title)
figure(3)
constellation_title = 'QPSK: With Noise';  
display_constellation_QPSK(r,constellation_title)
