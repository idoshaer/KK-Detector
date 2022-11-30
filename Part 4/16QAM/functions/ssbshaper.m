function [ssbshape,carrier] = ssbshaper(signal,beta,SR,Fs)
%%this function converts the signal to SSB configuration with carrier when
%%CSPR=10

B = (1+beta)*SR;
Ps = (abs(signal)).^2;
Ps_avg = mean(Ps);
Pc = 10*Ps_avg;
t = 0:(1/(Fs)):((length(signal)-1)/(Fs));
carrier = sqrt(Pc)*exp(-1j*pi*B*t);
ssbshape = signal + carrier;

end

