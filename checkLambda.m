a = 0.95;
b = 0.5;
nx = 5000;
x = sin(pi/2*linspace(-1,1,nx+2)); x(1) = []; x(end)=[];

gam0 = 3; gam1 = 4;
Lambda = -2*a./(1+a-b).*(gam0 + 2*(1-b)/(2+a-b)*gam1);

gam = @(xVar) gam0*weight(xVar,a,-b) + gam1*weight(xVar,a,1-b) + Lambda*weight(xVar,a-1,-b);

plot(x,gam(x));

lift = quadgk(gam,-1,1)
ylim([-5,5])