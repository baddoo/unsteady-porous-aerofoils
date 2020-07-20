function circFun = pFun2(xVar,gamCoefs,teZero,leSing,k)

jPInt = myJacobiP(numel(xVar),numel(gamCoefs)-3,1+teZero,2-leSing,xVar);

omega0 = exp(1i*k)*(2^(1+teZero-leSing)*beta(1-leSing,1+teZero))/(1+2i*k*beta(2-leSing,1));
omega1 = exp(1i*k)*(2^(2+teZero-leSing)*beta(2-leSing,1+teZero))/(1+2i*k*beta(2-leSing,1));

Lambda = -1i*k*2^(leSing-1)*exp(-1i*k)*(gamCoefs(1)*omega0 + gamCoefs(2)*omega1);

circFun =  gamCoefs(1)*2^(1+teZero-leSing)*myBeta((1+xVar)/2,1-leSing,1+teZero) ...
    + gamCoefs(2).*2^(2+teZero-leSing)*myBeta((1+xVar)/2,2-leSing,1+teZero) ...
    + sum(gamCoefs(3:end).'.*weight(xVar,1+teZero,2-leSing)./(-2*(1:(numel(gamCoefs)-2))).*jPInt,2) ...
    + Lambda*weight(xVar,0,2-leSing)/(2-leSing);

end
