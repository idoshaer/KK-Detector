function [inphase_1,inphase_2,inphase_3,quadrature_1,quadrature_2,quadrature_3] = decision64qam(noisySignal)
%%% this function uses the decision law for 64QAM symbols to estimate the
%%% correct bits from a noisy signal

%%% decision for the inphase bits
inphase_1(real(noisySignal)<=-6) = 0;
inphase_2(real(noisySignal)<=-6) = 0;
inphase_3(real(noisySignal)<=-6) = 0;

inphase_1(and(real(noisySignal) > -6,real(noisySignal) <= -4)) = 0;
inphase_2(and(real(noisySignal) > -6,real(noisySignal) <= -4)) = 0;
inphase_3(and(real(noisySignal) > -6,real(noisySignal) <= -4)) = 1;

inphase_1(and(real(noisySignal) > -4,real(noisySignal) <= -2)) = 0;
inphase_2(and(real(noisySignal) > -4,real(noisySignal) <= -2)) = 1;
inphase_3(and(real(noisySignal) > -4,real(noisySignal) <= -2)) = 1;

inphase_1(and(real(noisySignal) > -2,real(noisySignal) <= 0)) = 0;
inphase_2(and(real(noisySignal) > -2,real(noisySignal) <= 0)) = 1;
inphase_3(and(real(noisySignal) > -2,real(noisySignal) <= 0)) = 0;

inphase_1(and(real(noisySignal) > 0,real(noisySignal) <= 2)) = 1;
inphase_2(and(real(noisySignal) > 0,real(noisySignal) <= 2)) = 1;
inphase_3(and(real(noisySignal) > 0,real(noisySignal) <= 2)) = 0;

inphase_1(and(real(noisySignal) > 2,real(noisySignal) <= 4)) = 1;
inphase_2(and(real(noisySignal) > 2,real(noisySignal) <= 4)) = 1;
inphase_3(and(real(noisySignal) > 2,real(noisySignal) <= 4)) = 1;

inphase_1(and(real(noisySignal) > 4,real(noisySignal) <= 6)) = 1;
inphase_2(and(real(noisySignal) > 4,real(noisySignal) <= 6)) = 0;
inphase_3(and(real(noisySignal) > 4,real(noisySignal) <= 6)) = 1;

inphase_1(real(noisySignal) > 6) = 1;
inphase_2(real(noisySignal) > 6) = 0;
inphase_3(real(noisySignal) > 6) = 0;

%%% decision for the quadrature bits
quadrature_1(imag(noisySignal)<=-6) = 0;
quadrature_2(imag(noisySignal)<=-6) = 0;
quadrature_3(imag(noisySignal)<=-6) = 0;

quadrature_1(and(imag(noisySignal) > -6,imag(noisySignal) <= --4)) = 0;
quadrature_2(and(imag(noisySignal) > -6,imag(noisySignal) <= -4)) = 0;
quadrature_3(and(imag(noisySignal) > -6,imag(noisySignal) <= -4)) = 1;

quadrature_1(and(imag(noisySignal) > -4,imag(noisySignal) <= -2)) = 0;
quadrature_2(and(imag(noisySignal) > -4,imag(noisySignal) <= -2)) = 1;
quadrature_3(and(imag(noisySignal) > -4,imag(noisySignal) <= -2)) = 1;

quadrature_1(and(imag(noisySignal) > -2,imag(noisySignal) <= 0)) = 0;
quadrature_2(and(imag(noisySignal) > -2,imag(noisySignal) <= 0)) = 1;
quadrature_3(and(imag(noisySignal) > -2,imag(noisySignal) <= 0)) = 0;

quadrature_1(and(imag(noisySignal) > 0,imag(noisySignal) <= 2)) = 1;
quadrature_2(and(imag(noisySignal) > 0,imag(noisySignal) <= 2)) = 1;
quadrature_3(and(imag(noisySignal) > 0,imag(noisySignal) <= 2)) = 0;

quadrature_1(and(imag(noisySignal) > 2,imag(noisySignal) <= 4)) = 1;
quadrature_2(and(imag(noisySignal) > 2,imag(noisySignal) <= 4)) = 1;
quadrature_3(and(imag(noisySignal) > 2,imag(noisySignal) <= 4)) = 1;

quadrature_1(and(imag(noisySignal) > 4,imag(noisySignal) <= 6)) = 1;
quadrature_2(and(imag(noisySignal) > 4,imag(noisySignal) <= 6)) = 0;
quadrature_3(and(imag(noisySignal) > 4,imag(noisySignal) <= 6)) = 1;

quadrature_1(imag(noisySignal) > 6) = 1;
quadrature_2(imag(noisySignal) > 6) = 0;
quadrature_3(imag(noisySignal) > 6) = 0;
end

