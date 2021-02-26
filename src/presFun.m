function [fullP,ncP,qsP] = presFun(xVar,solStruct)

k = solStruct.k;

[fullV,ncV,qsV] = vortFun(xVar,solStruct);
[fullVint,ncVint] = vortInt(xVar,solStruct);

if k<inf
    
fullP = -2*(fullV + 1i*k*fullVint);

qsP = -2*qsV;

ncP = -2*(ncV + 1i*k*ncVint);

else
    
ncP = -2i*ncVint;
fullP = nan(size(ncP));
qsP = nan(size(ncP));

end
end