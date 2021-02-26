function w = weight(xVar,aVar,bVar)

w= (1-xVar).^aVar.*(1+xVar).^bVar;

end