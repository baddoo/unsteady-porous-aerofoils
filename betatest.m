xlin = linspace(-1,1); xlin(1)=[]; xlin(end)=[];
a = .15; b = .75;
z = (1-xlin)/2;
f1 = betainc(z,b,a)*beta(b,a);
f2 = (z.^b.*(1-z).^a + (a+b)*beta(b+1,a)*betainc(z,b+1,a))/b;
f2 = (z.^b.*(1-z).^a + (a+b)*(beta(b+1,a)-beta(a,b+1)*betainc(1-z,a,b+1)))/b;
f2 = (z.^b.*(1-z).^a + (a+b)*(beta(b+1,a)- ...
                              ((1-z).^(a).*(z).^(b+1) + (a+b+1)*beta(a+1,b+1)*betainc(1-z,a+1,b+1))/(a)))/b;
f2 = ((z).^(a).*(1-z).^(b) + (a+b)*beta(a+1,b)*betainc(z,a+1,b))/(a);

f1 = hypergeom([-a,1+b],1-a,.5*(1-xlin));
f2 = -a*2^(-a)*(1-xlin).^a.*myBeta(.5*(1-xlin),-a,-b);

plot(xlin,f1)
hold on
plot(xlin,f2)
hold off