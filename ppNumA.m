function pA = ppNumA(xVar,a,pCoefsA,psi0)

gam = 1/2 - acot(psi0)/pi;
tau2x = @(xVar)  1 + 2*(xVar-1)/(1-a);
delta = acot(psi0)/pi;

pA = sum(pCoefsA.'.*weight(tau2x(xVar),delta,gam).*myJacobiP(numel(xVar),numel(pCoefsA)-1,delta,gam,tau2x(xVar)),2);

end