function circFun = cFun2(xVar,gamCoefs,teZero,leSing,k)

jP = myJacobiP(numel(xVar),numel(gamCoefs)-2,teZero,1-leSing,xVar);

omega0 = exp(1i*k)*(2^(1+teZero-leSing)*beta(1-leSing,1+teZero))/(1+2i*k*beta(2-leSing,1));
omega1 = exp(1i*k)*(2^(2+teZero-leSing)*beta(2-leSing,1+teZero))/(1+2i*k*beta(2-leSing,1));

Lambda = -1i*k*2^(leSing-1)*exp(-1i*k)*(gamCoefs(1)*omega0 + gamCoefs(2)*omega1);

circFun =  gamCoefs(1)*(weight(xVar,teZero,-leSing)) ...
    + sum(gamCoefs(2:end).'.*weight(xVar,teZero,1-leSing).*jP,2) ...
    + Lambda*weight(xVar,0,1-leSing);

end