% Gives the right-hand side of the Jacobi polynomials orthogonality
% relation.

function f = jacobiOrth(n,a,b)

if n ==0
   
    f = 2^(a+b+1).*double(gamma(sym(a+1))).*double(gamma(sym(b+1)))./double(gamma(sym(a+b+2)));
    
else
    
    f =  2^(a+b+1)./(2*n+a+b+1).*double(gamma(sym(n+a+1))).*double(gamma(sym(n+b+1)))./double(gamma(sym(n+a+b+1)))./factorial(n);

end

end