% Applies recurrence formula to make MATLAB's incomplete Beta function valid for
% negative arguments

function B = myBeta(z,b,a)

if abs(imag(b))>0 | abs(imag(a))>0 | abs(imag(z))>0
    
     B = z.^b.*(1-z).^a/b.*hypergeom([a+b,1],b+1,z);
     
elseif  b<0 || -1<b || a<0  || -1<a

    B = (z.^b.*(1-z).^a + (a+b)*(beta(b+1,a+1)*(a+b+1)/a- ...
                              ((1-z).^(a).*(z).^(b+1) + (a+b+1)*beta(a+1,b+1)*betainc(1-z,a+1,b+1))/(a)))/b;
else

    B = betainc(z,b,a)*beta(b,a);    

end

end