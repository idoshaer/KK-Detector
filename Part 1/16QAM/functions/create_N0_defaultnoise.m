function [N0,default_noise,gamma_b,gamma_dB] = create_N0_defaultnoise(s,M,nos,sps,SNR_limit_value)
%%this function calculates the N0 for all SNR values and creates the default noise

%%calculation of the average energy per symbol and bit
epsilon_s = sum((abs(s)).^2,'all')/nos;
epsilon_b = epsilon_s/log2(M);

gamma_dB = 1:1:SNR_limit_value; %% defining the SNR values in dB
gamma_b = 10.^((gamma_dB)/10);  %% calculation of the SNR in linear values
N0 = epsilon_b./gamma_b;  %% calculation of the N0 for all SNR values
default_noise = randn(1,nos*sps)+1j*(randn(1,nos*sps)); %% creation of the default noise 

end

