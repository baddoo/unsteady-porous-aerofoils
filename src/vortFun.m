function [fullV,ncV,qsV] = vortFun(xVar,solStruct)

fullCoefs = solStruct.fullCoefs;
qsCoefs = solStruct.qsCoefs;
ncCoefs = solStruct.ncCoefs;
teZero = solStruct.teZero;
leSing = solStruct.leSing;
teZeroQS = solStruct.teZeroQS;
leSingQS = solStruct.leSingQS;
Gamma = solStruct.Gamma;
Theta = solStruct.Theta;
k = solStruct.k;

jP   = myJacobiP(numel(xVar),numel(fullCoefs)-2,teZero,1-leSing,xVar);

fullV = fullCoefs(1) * (weight(xVar,teZero,-leSing)) ...
    + sum(fullCoefs(2:end).' .* weight(xVar,teZero,1-leSing).*jP,2) ...
    - 1i*k*2^(leSing-1) * Gamma * weight(xVar,0,1-leSing);

ncV = ncCoefs(1) * (weight(xVar,teZero,-leSing)) ...
    + sum(ncCoefs(2:end).' .* weight(xVar,teZero,1-leSing).*jP,2) ...
    + Theta * weight(xVar,teZero-1,1-leSing);

if ~isnan(solStruct.teZeroQS)
jPQS = myJacobiP(numel(xVar),numel(fullCoefs)-2,teZeroQS,1-leSingQS,xVar);

qsV = qsCoefs(1) * (weight(xVar,teZeroQS,-leSingQS)) ...
    + sum(qsCoefs(2:end).' .* weight(xVar,teZeroQS,1-leSingQS).*jPQS,2);
else 
    qsV = nan;
end

end