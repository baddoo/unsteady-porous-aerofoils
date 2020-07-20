function pF = ppNumFCheb1(xVar,a,pCoefsF)

gam = 1e-10;
tau1x = @(xVar) -1 + 2*(xVar+1)/(1+a);

pF =  (pCoefsF(1)*weight(tau1x(xVar),gam,-1/2) ...
     + sum(pCoefsF(2:end).'.*weight(tau1x(xVar),gam,gam).*myJacobiP(numel(xVar),numel(pCoefsF)-2,gam,gam,tau1x(xVar)),2));
 
end