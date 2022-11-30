function [shifted_frequency] = KK_algorithm(downsampled_signal,Fs,beta,SR)
%%this function implements the KK algorithm and shows the spectrum graphs
%%of different stages in the process


sqrt_of_photocurrent = sqrt(downsampled_signal);
ln_of_sqrt_of_photocurrent = log(sqrt_of_photocurrent);
hilb = imag(hilbert(ln_of_sqrt_of_photocurrent));
exp_jphi = exp(1j*hilb);
downconverted_optical_field = sqrt_of_photocurrent.*exp_jphi;
figure(4)
title4 = 'downconverted optical field';
fft_plot_log(downconverted_optical_field,Fs/2,title4);

removed_carrier = downconverted_optical_field - mean(downconverted_optical_field);
t2 = 0:(1/(Fs/2)):((length(removed_carrier)-1)/(Fs/2));
B = (1+beta)*SR;
shifted_frequency = removed_carrier.*exp(-1j*pi*B*t2);
figure(5)
title5 = 'shifted frequency';
fft_plot_log(shifted_frequency,Fs/2,title5);

end

