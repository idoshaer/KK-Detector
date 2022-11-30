function [inphase_1,inphase_2,quadrature_1,quadrature_2] = decision16qam(noisySignal)
%%% this function uses the decision law for 16QAM symbols to estimate the
%%% correct bits from a noisy signal

%%% decision for the inphase bits
inphase_1(real(noisySignal)<=-2) = 0;
inphase_2(real(noisySignal)<=-2) = 0;

inphase_1(and(real(noisySignal) > -2,real(noisySignal) <= 0)) = 0;
inphase_2(and(real(noisySignal) > -2,real(noisySignal) <= 0)) = 1;

inphase_1(and(real(noisySignal) > 0,real(noisySignal) <= 2)) = 1;
inphase_2(and(real(noisySignal) > 0,real(noisySignal) <= 2)) = 1;

inphase_1(real(noisySignal) > 2) = 1;
inphase_2(real(noisySignal) > 2) = 0;

%%% decision for the quadrature bits
quadrature_1(imag(noisySignal)<=-2) = 0;
quadrature_2(imag(noisySignal)<=-2) = 0;

quadrature_1(and(imag(noisySignal) > -2,imag(noisySignal) <= 0)) = 0;
quadrature_2(and(imag(noisySignal) > -2,imag(noisySignal) <= 0)) = 1;

quadrature_1(and(imag(noisySignal) > 0,imag(noisySignal) <= 2)) = 1;
quadrature_2(and(imag(noisySignal) > 0,imag(noisySignal) <= 2)) = 1;

quadrature_1(imag(noisySignal) > 2) = 1;
quadrature_2(imag(noisySignal) > 2) = 0;

end

