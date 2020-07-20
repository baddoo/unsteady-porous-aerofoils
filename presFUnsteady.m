function pF = presFUnsteady(xVar,solStruct)

a = solStruct.junction;
midZero = solStruct.midZero;
tau1x = @(xVar) -1 + 2*(xVar+1)/(1+a);
nf = solStruct.nf;
forCoefs = solStruct.coefs(1:nf);
k = solStruct.k;

% PI coefficients
piCoefs = definePiCoefficients(solStruct);
Pi = piCoefs(1)*forCoefs(1) + piCoefs(2)*forCoefs(2);

vF =  forCoefs(1)*weight(tau1x(xVar),midZero,-1/2) ...
      +Pi/sqrt(2)*weight(tau1x(xVar),0,1/2) ...
     + sum(forCoefs(2:end).'.*weight(tau1x(xVar),midZero,1/2).*myJacobiP(numel(xVar),numel(forCoefs)-2,midZero,1/2,tau1x(xVar)),2);

vFInt=  (1+a)/2*(forCoefs(1).*myJacobiI(0,midZero,-1/2,tau1x(xVar)))...
       +(1+a)/2*Pi/sqrt(2).*myJacobiI(0,0,1/2,tau1x(xVar))...
       +(1+a)/2*sum(forCoefs(2:end).'.*myJacobiI(numel(forCoefs)-2,midZero,1/2,tau1x(xVar)),2);

pF = -2*(vF+ 1i*k*vFInt);
            
end