function [y] = fiber(s,Fv,D,L)

s_fft = fftshift(fft(s));

c = 3e17; %% [nm/s]
w = 2*pi*Fv;
lambda = 1550; %% [nm]
beta2 = -(((lambda)^2)*D)/(2*pi*c);

H_CD = exp((1j*beta2)/2*(w.^2)*L);

y_fft = s_fft.*H_CD;
y = ifft(ifftshift(y_fft));

end
