function  U = chebU(xVar,nVar)

U = sin((nVar+1).*acos(xVar))./sin(acos(xVar));

end
