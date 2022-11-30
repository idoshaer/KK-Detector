function [pulse_shaped_signal,srrc] = srrc_conv(signal,sps,span,beta)
%%this function convolutes the signal with a srrc shaped filter
%%and returns a pulse shaped signal without the tail elements

signal_length = length(signal);
srrc = rcosdesign(beta,span,sps,'sqrt');
pulse_shaped_signal_with_tails = conv(signal,srrc);
tail_elements = (span*sps)/2;
pulse_shaped_signal = pulse_shaped_signal_with_tails((tail_elements + 1):((signal_length) + tail_elements));


end

