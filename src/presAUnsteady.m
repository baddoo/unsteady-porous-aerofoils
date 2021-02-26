function pA = presAUnsteady(xVar,solStruct)

a = solStruct.junction;
midZero = solStruct.midZero;
teZero=solStruct.teZero;
tau2x = @(xVar)  1 + 2*(xVar-1)/(1-a);
nf = solStruct.nf;
aftCoefs = solStruct.coefs((nf+1):end);
forCoefs = solStruct.coefs(1:nf);
k = solStruct.k;

% PI coefficients
piCoefs = definePiCoefficients(solStruct);
Pi = piCoefs(1)*forCoefs(1) + piCoefs(2)*forCoefs(2);

lamCoefs = defineLambdaCoefs(solStruct);
finalLambda = (lamCoefs(1).*solStruct.coefs(1)...
             + lamCoefs(2).*solStruct.coefs(2)...
             + lamCoefs(nf+1).*aftCoefs(1));

vA = sum(aftCoefs.'.*weight(tau2x(xVar),teZero,midZero).*myJacobiP(numel(xVar),numel(aftCoefs)-1,teZero,midZero,tau2x(xVar)),2) ...
    + finalLambda*weight(tau2x(xVar),0,midZero) ...
    + Pi/2^teZero*weight(tau2x(xVar),teZero,0);

aftInt = forCoefs(1)*2^(.5+midZero)*double(beta(.5,sym(1+midZero))) ...
        +forCoefs(2)*2^(1.5+midZero)*double(beta(1.5,sym(1+midZero))) ...
        +Pi/sqrt(2)*2^(3/2)*beta(1,3/2);

vAInt = (1+a)/2*aftInt + (1-a)/2*sum(aftCoefs.'.*myJacobiI(numel(aftCoefs)-1,teZero,midZero,tau2x(xVar)),2) ...
                       + (1-a)/2*finalLambda*myJacobiI(0,0,midZero,tau2x(xVar))...
        +(1-a)/2*Pi/2^teZero*myJacobiI(0,teZero,0,tau2x(xVar));

pA = -2*(vA + 1i*k*vAInt);
                               
end