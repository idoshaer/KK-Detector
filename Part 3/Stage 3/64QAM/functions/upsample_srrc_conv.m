function [pulse_shaped_signal,srrc,tail_elements] = upsample_srrc_conv(s,nos,sps,span,beta)
%%this function upsamples the given signal and convolutes it with a srrc shaped filter
%%the function returns a pulse shaped signal without the tail elements

tail_elements = (span*sps)/2;
s_upsampled = upsample(s,sps);
srrc = rcosdesign(beta,span,sps,'sqrt');
pulse_shaped_signal_with_tails = conv(s_upsampled,srrc);
pulse_shaped_signal = pulse_shaped_signal_with_tails((tail_elements + 1):((nos*sps) + tail_elements));

end

