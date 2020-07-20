function pF = ppNumF(xVar,a,pCoefsF,psi0)

gam = 1/2 - acot(psi0)/pi;
tau1x = @(xVar) -1 + 2*(xVar+1)/(1+a);

pF =  (pCoefsF(1)*weight(tau1x(xVar),gam,-1/2) ...
     + sum(pCoefsF(2:end).'.*weight(tau1x(xVar),gam,1/2).*myJacobiP(numel(xVar),numel(pCoefsF)-2,gam,1/2,tau1x(xVar)),2));
 
end