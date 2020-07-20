function piCoefs = definePiCoefficients(solStruct)

a = solStruct.junction;
k = solStruct.k;
midZero = solStruct.midZero;
piCoefs = zeros(1,solStruct.nf + solStruct.na);

piPreFac = -1i*k*(1+a)/(2*(1+1i*k*(1+a)*beta(3/2,1)));
piCoefs(1) = 2^(.5 +midZero)*double(beta(.5,sym(1+midZero)))*piPreFac;
piCoefs(2) = 2^(1.5+midZero)*double(beta(1.5,sym(1+midZero)))*piPreFac;

end