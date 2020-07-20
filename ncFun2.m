function pF = ncFun(xVar,gamCoefs,alpha,beta)

jP = myJacobiP(numel(xVar),numel(gamCoefs)-2,alpha,1-beta,xVar);

pF =  gamCoefs(1)*(weight(xVar,alpha,-beta) - 2*alpha/(1+alpha-beta).*weight(xVar,alpha-1,-beta)) ...
    + gamCoefs(2)*(weight(xVar,alpha,1-beta)- 4*alpha*(1-beta)/(1+alpha-beta)/(2+alpha-beta).*weight(xVar,alpha-1,-beta)) ...
    + sum(gamCoefs(3:end).'.*weight(xVar,alpha,1-beta).*jP(:,2:end),2);

end