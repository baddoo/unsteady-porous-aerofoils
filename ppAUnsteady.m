function pA = ppAUnsteady(xVar,solStruct)

a = solStruct.junction;
midZero = solStruct.midZero;
teZero=solStruct.teZero;
tau2x = @(xVar)  1 + 2*(xVar-1)/(1-a);
aftCoefs = solStruct.coefs((solStruct.nf+1):end);
k = solStruct.k;

preFac = -1i*k*2^(-midZero)/(1 + 1i*k*(1-a)*beta(1,1+midZero));
gam0iWC = preFac*(1+a)/2*2^(.5+midZero)*beta(.5,1+midZero)            ;
gam1iWC = preFac*(1+a)/2*2^(1.5+midZero)*beta(1.5,1+midZero)          ;
gam0pWC = preFac*(1-a)/2*2^(1+teZero+midZero)*beta(1+teZero,1+midZero);

pA = sum(aftCoefs.'.*weight(tau2x(xVar),teZero,midZero).*myJacobiP(numel(xVar),numel(aftCoefs)-1,teZero,midZero,tau2x(xVar)),2) ...
    +weight(tau2x(xVar),0,midZero).*(gam0iWC.*solStruct.coefs(1)...
                               + gam1iWC.*solStruct.coefs(2)...
                               + gam0pWC.*aftCoefs(1));

end