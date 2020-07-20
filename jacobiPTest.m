x = linspace(-1,1);
a= -.5; b= -.55;
n = 15;
exact = jacobiP(n,a,b,x);
numPols = myJacobiP(numel(x),n+1,a,b,x.').';
approx = numPols(n+1,:);

plot(x,abs(exact-approx));