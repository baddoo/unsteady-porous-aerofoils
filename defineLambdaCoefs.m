function lamCoefs = defineLambdaCoefs(solStruct)

a = solStruct.junction;
k = solStruct.k;
midZero = solStruct.midZero;
teZero = solStruct.teZero;
nf = solStruct.nf; na = solStruct.na;
piCoefs = definePiCoefs(solStruct);
lamCoefs = zeros(1,na+nf);

preFac = -1i*k*2^(-midZero)/(2*(1 + 1i*k*(1-a)*double(beta(1,sym(1+midZero)))));
lamCoefs(1) = preFac*((1+a)*(2*piCoefs(1)*beta(1,1.5) ...
                      +2^(.5+midZero)*double(beta(.5,sym(1+midZero))))...
               +(1-a)*2*piCoefs(1)*double(beta(1,sym(1+teZero))));
lamCoefs(2) = preFac*((1+a)*(2*piCoefs(2)*beta(1,1.5) ...
                      +2^(1.5+midZero)*double(beta(1.5,sym(1+midZero))))...
               +(1-a)*2*piCoefs(2)*double(beta(1,sym(1+teZero))));
lamCoefs(nf+1) = preFac*(1-a)*2^(1+teZero+midZero)*double(beta(sym(1+teZero),sym(1+midZero)));

end