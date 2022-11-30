function [] = BER_delta_L_SR_plot(delta_L,BER_matrix,mod_title)
%%% this function plots the BER for different length deviations of
%%% the optical fiber in 3 different symbol rates.

subplot(2,2,1)
semilogy(delta_L,BER_matrix(1,:));
title('SR = 15 [Gbaud]');
xlabel('Delta L [m]');
ylabel('BER');
grid on
subplot(2,2,2)
semilogy(delta_L,BER_matrix(2,:));
title('SR = 30 [Gbaud]');
xlabel('Delta L [m]');
ylabel('BER');
grid on
subplot(2,2,3)
semilogy(delta_L,BER_matrix(3,:));
title('SR = 60 [Gbaud]');
xlabel('Delta L [m]');
ylabel('BER');
grid on
subplot(2,2,4)
semilogy(delta_L,BER_matrix(1,:),'b',delta_L,BER_matrix(2,:),'g',delta_L,BER_matrix(3,:),'r');
title('Comparison Between Different Symbol Rates');
xlabel('Delta L [m]');
ylabel('BER');
grid on
legend('SR = 15 [Gbaud]','SR = 30 [Gbaud]','SR = 60 [Gbaud]');

sgtitle(mod_title);

end

