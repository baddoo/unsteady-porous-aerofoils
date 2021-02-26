% The definite integral of the vorticity distribution from -1 to xVar

function [fullVint,ncVint,qsVint] = vortInt(xVar,solStruct)

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

jPInt = myJacobiP(numel(xVar),numel(fullCoefs)-3,1+teZero,2-leSing,xVar);
jPIntQS = myJacobiP(numel(xVar),numel(fullCoefs)-3,1+teZeroQS,2-leSingQS,xVar);

fullVint =  fullCoefs(1)*2^(1+teZero-leSing)*myBeta((1+xVar)/2,1-leSing,1+teZero) ...
         + fullCoefs(2).*2^(2+teZero-leSing)*myBeta((1+xVar)/2,2-leSing,1+teZero) ...
         + sum(fullCoefs(3:end).'.*weight(xVar,1+teZero,2-leSing)./(-2*(1:(numel(fullCoefs)-2))).*jPInt,2) ...
         - 1i*k*Gamma*2^(leSing-1)*weight(xVar,0,2-leSing)/(2-leSing);

ncVint =   ncCoefs(1)*2^(1+teZero-leSing)*myBeta((1+xVar)/2,1-leSing,1+teZero) ...
         + ncCoefs(2)*2^(2+teZero-leSing)*myBeta((1+xVar)/2,2-leSing,1+teZero) ...
         + sum(ncCoefs(3:end).'.*weight(xVar,1+teZero,2-leSing)./(-2*(1:(numel(ncCoefs)-2))).*jPInt,2) ...
         + Theta*2^(1+teZero-leSing)*myBeta((1+xVar)/2,2-leSing,teZero);     

 if ~isnan(solStruct.teZeroQS)  
qsVint =  qsCoefs(1)*2^(1+teZeroQS-leSingQS)*myBeta((1+xVar)/2,1-leSingQS,1+teZeroQS) ...
     + qsCoefs(2).*2^(2+teZeroQS-leSingQS)*myBeta((1+xVar)/2,2-leSingQS,1+teZeroQS) ...
     + sum(qsCoefs(3:end).'.*weight(xVar,1+teZeroQS,2-leSingQS)./(-2*(1:(numel(qsCoefs)-2))).*jPIntQS,2);
else 
    qsVint = nan;
end

end
