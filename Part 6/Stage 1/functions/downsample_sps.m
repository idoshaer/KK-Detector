function [downsampled_signal,sps_new] = downsample_sps(signal,sps,decimation_index,offset)
%%this function downsamples the given signal and returns the downsampled signal and the new sps

downsampled_signal = downsample(signal,decimation_index,offset);
sps_new = sps/decimation_index;

end

