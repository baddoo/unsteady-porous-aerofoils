function I = myJacobiI(nVar,aVar,bVar,xVar)

xVar = reshape(xVar,[numel(xVar),1]);

    I(:,1) = 2^aVar.*(1+xVar).^(1+bVar)./(1+bVar).*hypergeom([-aVar,1+bVar],2+bVar,(1+xVar)/2);
    %I(:,1) = 2^(1+aVar+bVar).*betainc((1+xVar)/2,1+bVar,1+aVar).*beta(1+bVar,1+aVar);
    
if nVar>0

    I(:,2:(nVar+1)) = -weight(xVar,aVar+1,bVar+1)./(2.*(1:(nVar))).*myJacobiP( numel(xVar), nVar-1, aVar+1, bVar+1, xVar );
    
end

end