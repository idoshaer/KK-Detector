function [down_sampled] = srrc_conv_downsample(r,srrc,tail_elements,signal_length,sps,offset)
%%this function convolutes the given signal with a srrc shaped filter and downsamples it
%%the function returns a noisy pulse shaped signal without the tail elements

pulse_shaped_signal_noisy_with_tails = conv(r,srrc);
pulse_shaped_signal_noisy = pulse_shaped_signal_noisy_with_tails((tail_elements + 1):((signal_length) + tail_elements));
down_sampled = downsample(pulse_shaped_signal_noisy,sps,offset);

end
