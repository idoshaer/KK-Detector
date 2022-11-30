function display_constellation_64QAM(IQ,constellation_title)

%for display of QPSK constellation diagram with color density
%Ix=real part of teh symbol vector
%Qx=imaginary aprt of the symbol vector

Ix=real(IQ);
Qx=imag(IQ);
data=[Ix',Qx'];
count=hist2d(data,-9.5:0.05:9.5,-9.5:0.05:9.5);
imagesc(-9.5:0.05:9.5,-9.5:0.05:9.5,count);
set(gca,'YDir','normal')
colormap([1,1,1;jet(128)]) 
title([constellation_title,' Constellation'])
xlabel('In-phase')
ylabel('Quadrature')
grid on