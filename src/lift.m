% Calculates the lift using the orthogonality relations. Integrating by
% parts allows the lift to be expressed in terms of the circulation and
% another integral.

function [fullLift,ncLift,qsLift] = lift(solStruct)

Theta = solStruct.Theta;
Gamma = solStruct.Gamma;
fullCoefs = solStruct.fullCoefs;
ncCoefs = solStruct.ncCoefs;
teZero = solStruct.teZero;
leSing = solStruct.leSing;
k = solStruct.k;

% Jacobi orthogonality constants
jOXa = jacobiOrthX(0,teZero,-leSing);
jOXb = jacobiOrthX([0,1],teZero,1-leSing);

% Calculate the circulation
[fullCirc,ncCirc,qsCirc] = circulation(solStruct);
 
% Calculate the second integrals
fullInt =  fullCoefs(1) * jOXa + sum(fullCoefs(2:3).' .* jOXb) ...
      - 1i*k*2^(leSing-1) * Gamma * jacobiOrthX(0,0,1-leSing);
    
ncInt =  ncCoefs(1) * jOXa + sum(ncCoefs(2:3).' .* jOXb)...
      + Theta * jacobiOrthX(0,teZero-1,1-leSing);
  
% Combine  
fullLift = -2 * ( fullCirc + 1i*k * (fullCirc - fullInt) ) ; 
qsLift = -2 * (qsCirc) ; 
ncLift = -2 * ( ncCirc + 1i*k * (ncCirc - ncInt) ) ; 

end