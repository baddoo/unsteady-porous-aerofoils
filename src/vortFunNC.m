function circFun = vortFunNC(xVar,gamCoefs,teZero,leSing,k)

jP = myJacobiP(numel(xVar),numel(gamCoefs)-2,teZero,1-leSing,xVar);

denom = -2^(1+teZero-leSing)*double(beta(sym(teZero),sym(2-leSing)));
omega0 = (2^(1+teZero-leSing)*double(beta(sym(1-leSing),sym(1+teZero))))/denom;
omega1 = (2^(2+teZero-leSing)*double(beta(sym(2-leSing),sym(1+teZero))))/denom;    

omeg = omega0*gamCoefs(1) + omega1*gamCoefs(2);

circFun =  gamCoefs(1)*(weight(xVar,teZero,-leSing)) ...
    + sum(gamCoefs(2:end).'.*weight(xVar,teZero,1-leSing).*jP,2) ...
    + omeg*weight(xVar,teZero-1,1-leSing);

end