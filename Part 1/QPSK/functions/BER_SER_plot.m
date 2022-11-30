function [] = BER_SER_plot(mod_title,gamma_dB,BER,theoretical_BER,SER,theoretical_SER,SNR_limit_value,y_limits)
%%this function plots the graph of the theoretical and the numerical BER and SER as a
%%function of the SNR values in dB.

subplot(2,1,1)
BER_plot = semilogy(gamma_dB,BER,'g',gamma_dB,theoretical_BER,'--ob');
BER_plot(1).LineWidth = 2;
BER_plot(2).LineWidth = 2;
xlim([1 SNR_limit_value])
ylim(y_limits)
xlabel('SNR [dB]')
ylabel('BER')
title([mod_title,' BER'])
legend('Numeric BER','Theoretical BER')
grid on


subplot(2,1,2)
SER_plot = semilogy(gamma_dB,SER,'g',gamma_dB,theoretical_SER,'--ob');
SER_plot(1).LineWidth = 2;
SER_plot(2).LineWidth = 2;
xlim([1 SNR_limit_value])
ylim(y_limits)
xlabel('SNR [dB]')
ylabel('SER')
title([mod_title,' SER'])
legend('Numeric SER','Theoretical SER')
grid on

end

