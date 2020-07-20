% Checks that our Jacobi integral function is correct.

x = linspace(-1,1,1000); x(1)=[]; x(end)=[];
n = 30;
a = .25;
b = .125;

exact = myJacobiI(x,a,b,n);
jac = myJacobiP( numel(x), n, a, b, x.' ).';
jac3 = permute(jac,[3,2,1]);
approx3 = cumtrapz(x,weight(x,a,b).*jac3,2);
approx = permute(approx3,[3,2,1]);

semilogy(x,abs(exact-approx))