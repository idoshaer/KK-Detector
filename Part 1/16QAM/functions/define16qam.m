function [symbols,inphase_1,inphase_2,quadrature_1,quadrature_2] = define16qam(numberOfSignals)
%%this function creates bit streams to define gray coded 16QAM symbols
%%with separation of inphase and quadrature

%%creation of the inphase logical bit stream source randomly
inphase_1 = rand(1,numberOfSignals)>0.5;
inphase_2 = rand(1,numberOfSignals)>0.5;

%%creation of the quadrature logical bit stream source randomly
quadrature_1 = rand(1,numberOfSignals)>0.5;
quadrature_2 = rand(1,numberOfSignals)>0.5;

indexI_00 = and(inphase_1 == 0,inphase_2 == 0);
indexI_10 = and(inphase_1 == 1,inphase_2 == 0);
indexI_11 = and(inphase_1 == 1,inphase_2 == 1);
indexI_01 = and(inphase_1 == 0,inphase_2 == 1);

indexQ_00 = and(quadrature_1 == 0,quadrature_2 == 0);
indexQ_10 = and(quadrature_1 == 1,quadrature_2 == 0);
indexQ_11 = and(quadrature_1 == 1,quadrature_2 == 1);
indexQ_01 = and(quadrature_1 == 0,quadrature_2 == 1);

%%assigning numerical values to each of the inphase bit combinations
I(indexI_00) = -3;
I(indexI_01) = -1;
I(indexI_11) = 1;
I(indexI_10) = 3;

%%assigning numerical values to each of the quadrature bit combinations
Q(indexQ_00) = -3;
Q(indexQ_01) = -1;
Q(indexQ_11) = 1;
Q(indexQ_10) = 3;


%%creation of the 16QAM symbols with seperation of the inphase and the quadrature
symbols = I+1j*Q;

end

