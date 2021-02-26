function pF = ppFUnsteady(xVar,solStruct)

a = solStruct.junction;
midZero = solStruct.midZero;
tau1x = @(xVar) -1 + 2*(xVar+1)/(1+a);
forCoefs = solStruct.coefs(1:(solStruct.nf));

pF =  (forCoefs(1)*weight(tau1x(xVar),midZero,-1/2) ...
     + sum(forCoefs(2:end).'.*weight(tau1x(xVar),midZero,1/2).*myJacobiP(numel(xVar),numel(forCoefs)-2,midZero,1/2,tau1x(xVar)),2));
 
end