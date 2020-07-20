

z = linspace(-1,1); z(1) = []; z(end) = [];

b = .5+0.2i; a = 0.9+.1i;

n = 1;

jacobiP(n,a,b,z);

norm(jacobiP(n,a,b-1,z)-jacobiP(n,a-1,b,z) - jacobiP(n-1,a,b,z),'inf')

return
tic
 b1 = myBeta(z,b,a);
 toc
 %b2 = betainc(z,b,a)*beta(b,a);
 tic
 b2 = z.^b.*(1-z).^a/b.*hypergeom([a+b,1],b+1,z);
 toc
 error = norm(b1-b2,'inf')
 