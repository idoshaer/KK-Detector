function [symbols,inphase_1,inphase_2,inphase_3,quadrature_1,quadrature_2,quadrature_3] = define64qam(numberOfSignals)
%%this function creates bit streams to define gray coded 64QAM symbols
%%with separation of inphase and quadrature

%%creation of the inphase logical bit stream source randomly
inphase_1 = rand(1,numberOfSignals)>0.5;
inphase_2 = rand(1,numberOfSignals)>0.5;
inphase_3 = rand(1,numberOfSignals)>0.5;

%%creation of the quadrature logical bit stream source randomly
quadrature_1 = rand(1,numberOfSignals)>0.5;
quadrature_2 = rand(1,numberOfSignals)>0.5;
quadrature_3 = rand(1,numberOfSignals)>0.5;


indexI_000 = inphase_1 == 0 & inphase_2 == 0 & inphase_3 == 0;
indexI_001 = inphase_1 == 0 & inphase_2 == 0 & inphase_3 == 1;
indexI_011 = inphase_1 == 0 & inphase_2 == 1 & inphase_3 == 1;
indexI_010 = inphase_1 == 0 & inphase_2 == 1 & inphase_3 == 0;
indexI_110 = inphase_1 == 1 & inphase_2 == 1 & inphase_3 == 0;
indexI_111 = inphase_1 == 1 & inphase_2 == 1 & inphase_3 == 1;
indexI_101 = inphase_1 == 1 & inphase_2 == 0 & inphase_3 == 1;
indexI_100 = inphase_1 == 1 & inphase_2 == 0 & inphase_3 == 0;

indexQ_000 = quadrature_1 == 0 & quadrature_2 == 0 & quadrature_3 == 0;
indexQ_001 = quadrature_1 == 0 & quadrature_2 == 0 & quadrature_3 == 1;
indexQ_011 = quadrature_1 == 0 & quadrature_2 == 1 & quadrature_3 == 1;
indexQ_010 = quadrature_1 == 0 & quadrature_2 == 1 & quadrature_3 == 0;
indexQ_110 = quadrature_1 == 1 & quadrature_2 == 1 & quadrature_3 == 0;
indexQ_111 = quadrature_1 == 1 & quadrature_2 == 1 & quadrature_3 == 1;
indexQ_101 = quadrature_1 == 1 & quadrature_2 == 0 & quadrature_3 == 1;
indexQ_100 = quadrature_1 == 1 & quadrature_2 == 0 & quadrature_3 == 0;

%%assigning numerical values to each of the inphase bit combinations
I(indexI_000) = -7;
I(indexI_001) = -5;
I(indexI_011) = -3;
I(indexI_010) = -1;
I(indexI_110) = 1;
I(indexI_111) = 3;
I(indexI_101) = 5;
I(indexI_100) = 7;

%%assigning numerical values to each of the quadrature bit combinations
Q(indexQ_000) = -7;
Q(indexQ_001) = -5;
Q(indexQ_011) = -3;
Q(indexQ_010) = -1;
Q(indexQ_110) = 1;
Q(indexQ_111) = 3;
Q(indexQ_101) = 5;
Q(indexQ_100) = 7;

%%creation of the 64QAM symbols with seperation of the inphase and the quadrature
symbols = I+1j*Q;

end

