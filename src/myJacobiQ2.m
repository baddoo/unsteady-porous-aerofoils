function v = myJacobiQ2 (n, aVar, bVar, xVar )

m = size(xVar,1);

  if ( aVar <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'J_POLYNOMIAL - Fatal error!\n' );
    fprintf ( 1, '  Illegal input value of ALPHA = %f\n', aVar );
    fprintf ( 1, '  But ALPHA must be greater than -1.\n' );
    error ( 'J_POLYNOMIAL - Fatal error!' );
  end
 
  if ( bVar <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'J_POLYNOMIAL - Fatal error!\n' );
    fprintf ( 1, '  Illegal input value of BETA = %f\n', bVar );
    fprintf ( 1, '  But BETA must be greater than -1.\n' );
    error ( 'J_POLYNOMIAL - Fatal error!' );
  end
  
  if ( n < 0 )
    v = [];
    return
  end
 
  if all(-1<xVar) && all(xVar<1)

         v(1:m,1) = pi.*cot(pi*aVar)...
              -2.^(aVar).*double(beta(sym(bVar+1),sym(aVar+1)))*(aVar+bVar+1)/aVar./(1-xVar).^aVar.*-aVar.*2^(-aVar).*(1-xVar).^aVar.*myBeta(.5*(1-xVar),-aVar,-bVar);

  else
    
      prefac = -2^(aVar+bVar+1)*double(gamma(sym(aVar+1)))*double(gamma(sym(bVar+1)))/double(gamma(sym(aVar+bVar+2)))./weight(xVar,aVar,bVar);
      
        if all(abs(xVar-1)>2)  
      
            v(1:m,1) = prefac./(xVar-1).*hypergeom([aVar+1,1],aVar+bVar+2,2./(1-xVar));
            
         elseif all(abs(xVar+1)>2)      
        
            v(1:m,1) = prefac./(xVar+1).*hypergeom([bVar+1,1],aVar+bVar+2,2./(1+xVar));
        
        end
  end
  
  
  if ( n == 0 )
    return
  end
  
  % Recurrence relation:
  B0 = (aVar.^2-bVar.^2)./(aVar+bVar)./(aVar+bVar+2);
  A0 = 2*(aVar+bVar+1)./(aVar+bVar+1)./(aVar+bVar+2);
  M0 = 2^(aVar+bVar+1).*double(gamma(sym(aVar+1))).*double(gamma(sym(bVar+1)))./double(gamma(sym(aVar+bVar+2)));
  
  v(1:m,2) = (M0+(xVar+B0).*v(1:m,1).*weight(xVar,aVar,bVar))/A0./weight(xVar,aVar,bVar);

  for i = 2 : n

    c1 = 2 * i * ( i + aVar + bVar ) * ( 2 * i - 2 + aVar + bVar );

    c2 = ( 2 * i - 1 + aVar + bVar ) * ( 2 * i + aVar + bVar ) ...
      * ( 2 * i - 2 + aVar + bVar );

    c3 = ( 2 * i - 1 + aVar + bVar ) * ( aVar + bVar ) * ( aVar - bVar );

    c4 = - 2 * ( i - 1 + aVar ) * ( i - 1 + bVar )  * ( 2 * i + aVar + bVar );

    v(1:m,i+1) = ( ( c3 + c2 * xVar(1:m) ) .* v(1:m,i) + c4 * v(1:m,i-1) ) / c1;

  end

end