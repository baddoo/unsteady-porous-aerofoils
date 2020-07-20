function pF = cpNum(xVar,pCoefsF,alpha,beta)

pF =  pCoefsF(1)*weight(xVar,alpha,-beta) ...
     + sum(pCoefsF(2:end).'.*weight(xVar,alpha,1-beta).*myJacobiP(numel(xVar),numel(pCoefsF)-2,alpha,1-beta,xVar),2);

end