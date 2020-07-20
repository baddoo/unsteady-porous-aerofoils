function p = presFunNC(xVar,struct)

coefs = struct.coefs;
teZero = struct.teZero;
leSing = struct.leSing;
k = struct.k;

p = -2*(vortFunNC(xVar,coefs,teZero,leSing,k) ...
  + 1i*k*vortIntNC(xVar,coefs,teZero,leSing,k));

end