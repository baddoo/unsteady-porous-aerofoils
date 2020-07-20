function pA = ppNumACheb2(xVar,a,pCoefsA)

gam = 1/2;
tau2x = @(xVar)  1 + 2*(xVar-1)/(1-a);
delta = 1/2;

pA = sum(pCoefsA.'.*weight(tau2x(xVar),delta,gam).*myJacobiP(numel(xVar),numel(pCoefsA)-1,delta,gam,tau2x(xVar)),2);

end