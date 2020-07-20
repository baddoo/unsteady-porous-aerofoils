% Gives the right-hand side of the integral of a weighted Jacobi polynomial
% multiplied by X.

function f = jacobiOrthX(n,a,b)

loc0 = find(n==0);
loc1 = find(n==1);

f = zeros(size(n));
f(loc0) = (1-2*(a+1)/(a+b+2)) * jacobiOrth(0,a,b);
f(loc1) = 2/(a+b+2)*jacobiOrth(1,a,b);

end