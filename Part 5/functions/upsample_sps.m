function [upsampled_signal,sps_new] = upsample_sps(signal,sps,interpolation_index)
%%this function upsamples the given signal and returns the upsampled signal and the new sps

upsampled_signal = upsample(signal,interpolation_index);
sps_new = sps*interpolation_index;

end

