function circFun = vortIntNC(xVar,gamCoefs,teZero,leSing,k)

jPInt = myJacobiP(numel(xVar),numel(gamCoefs)-3,1+teZero,2-leSing,xVar);

denom = -2^(1+teZero-leSing)*double(beta(sym(teZero),sym(2-leSing)));
omega0 = (2^(1+teZero-leSing)*double(beta(sym(1-leSing),sym(1+teZero))))/denom;
omega1 = (2^(2+teZero-leSing)*double(beta(sym(2-leSing),sym(1+teZero))))/denom;    

omeg = omega0*gamCoefs(1) + omega1*gamCoefs(2);

circFun =  gamCoefs(1)*2^(1+teZero-leSing)*myBeta((1+xVar)/2,1-leSing,1+teZero) ...
    + gamCoefs(2).*2^(2+teZero-leSing)*myBeta((1+xVar)/2,2-leSing,1+teZero) ...
    + sum(gamCoefs(3:end).'.*weight(xVar,1+teZero,2-leSing)./(-2*(1:(numel(gamCoefs)-2))).*jPInt,2) ...
    + omeg*myJacobiI(0,teZero-1,1-leSing,xVar);

end
