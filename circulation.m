% Calculates the lift using the orthogonality relations. Integrating by
% parts allows the lift to be expressed in terms of the circulation and
% another integral.

function [fullCirc,ncCirc,qsCirc] = circulation(solStruct)

Gamma = solStruct.Gamma;
Theta = solStruct.Theta;
fullCoefs = solStruct.fullCoefs;
ncCoefs = solStruct.ncCoefs;
qsCoefs = solStruct.qsCoefs;

teZero = solStruct.teZero;
leSing = solStruct.leSing;
teZeroQS = solStruct.teZeroQS;
leSingQS = solStruct.leSingQS;

% Orthogonality constants for Jacobi polynomials
jOA = jacobiOrth(0,teZero,-leSing);
jOB = jacobiOrth(0,teZero,1-leSing);
jOAQS = jacobiOrth(0,teZeroQS,-leSingQS);
jOBQS = jacobiOrth(0,teZeroQS,1-leSingQS);

% Define common parts of circulation
ncC0 = ncCoefs(1)*jOA + ncCoefs(2).*jOB;
qsC0 = qsCoefs(1)*jOAQS + qsCoefs(2).*jOBQS;

% Define other parts of circualation
ncC1 = Theta*jacobiOrth(0,teZero-1,1-leSing);

% Combine
fullCirc = Gamma;
qsCirc = qsC0;
ncCirc = ncC0 + ncC1;
 
end