function [f B]=inbeta(z,a,b,N)
%F=INBETA(X,P,Q) In complete beta function for complex x
% there are branch cuts [-infinity,0] and [1 infinity] on the real axis
% f = int(y.^(a-1).*(1-y).^(b-1),0..x);
% B is the complete integral
% The calculation is based in Muir's continued fraction expansion
% Biometrika Vol. 22 284--297 (1930-31)
% discussed in
% L.A. Aroian Annals of Mathematical Statistics Vol. 12 218-223 (1941) 
% and correection vol 30 1265 (1959)
%
% for large $z$ we use a Taylor expansion in powers of 1/z 
%
% The number of terms N can be adjusted
%
% Copyright Jim McElwaine 2010
%
% Call without arguments to run a test with random a and b in the
% complex plane
% The maximum error in the derivative is returned as a test
if nargin<2
  test_inbeta;
  return;
end
if nargin<4
  N=[];
end
if isempty(N)
  N=20;
end
zmax = 5;
s0 = (a+1)./(a+b+2);
B = double(beta(sym(a),sym(b)));
az = abs(z);
f1 = find(real(z)<=s0 & az<zmax);
f2 = find(real(z)>s0  & az<zmax);
f3 = find(az>=zmax);
f = repmat(NaN,size(z));
ff = repmat(NaN,size(z));
f(f1) = inbeta3(z(f1),a,b,N);
f(f2) = B-inbeta3(1-z(f2),b,a,N);
f(f3) = inbeta2(z(f3),a,b,B,N);
return;
% Taylor series in z  
% expansion about 0
function f=inbeta1(z,a,b,N)
C=z.^a/a;
f=C;
for n=1:N
  C = z.*(C*(n-b)*(n+a-1)/(n*(n+a)));
  f = f+C;
  if max(abs(C))<tol
    break;
  end
end
% Taylor series in 1/z  
% expansion about infinity
function f=inbeta2(z,a,b,B,N)
sz = sign(imag(z));
C = -z.^(a+b-1)/(a+b-1).*exp(-i*b*pi*sz);
z = 1./z;
f = C;
for n=1:N
  C = z.*C*((n-b)*(n-a-b)/(n*(n+1-a-b)));
  f = f+C;
end
f=f+B*sin(pi*b)*exp(i*pi*a*sz)/sin(pi*a+pi*b);
% Continued fraction expansion about z=0
function f=inbeta3(z,a,b,N)
f=0;
for k=N:-1:1
  f = k*(b-k)*z./((a+2*k-1).*(a+2*k).*(1+f));
  j=k-1;
  f = -(a+j)*(a+b+j)*z./((a+2*j).*(a+2*j+1).*(1+f));
end
f=z.^a.*(1-z).^b./(a.*(1+f));
function test_inbeta
a=3*rand(1);
b=3*rand(1);
x=linspace(-10,10,101);
y=linspace(-10,10,100);
[yg xg]=ndgrid(y,x);
  
z = xg+i*yg;
f  = inbeta(z,a,b);
% Choose a random direction for differentiation
% imag(z) is away from the branch cuts so no problems
dx =  eps^(1/6)*exp(i*2*pi*rand(1)); % We're using a 6th order accurate differentiation scheme
		 
f6 = inbeta(z+3*dx,a,b);
f5 = inbeta(z+2*dx,a,b);
f4 = inbeta(z+1*dx,a,b);
f3 = inbeta(z-1*dx,a,b);
f2 = inbeta(z-2*dx,a,b);
f1 = inbeta(z-3*dx,a,b);
%df1 = (f4-f3)(2*dx);                  % 2nd order
%df1 = (8*(f4-f3)+(f2-f5))/(12*dx);    % 4th order
df1 = (45*(f4-f3)+9*(f2-f5)+(f6-f1))/(60*dx); % 3rd order
df2 = z.^(a-1).*(1-z).^(b-1);
disp(sprintf('max error %g',max(abs(df1(:)-df2(:)))));
nc = 30;
clf;
subplot(2,2,1);
contour(x,y,angle(f),nc);hold('on');plot([[x(1);0] [1;x(end)]],zeros(2),'k');
title(sprintf('arg, a = %5.3f, b = %5.3f',a,b));
subplot(2,2,2);
contour(x,y,abs(f),nc);hold('on');plot([[x(1);0] [1;x(end)]],zeros(2),'k');
title('abs');
subplot(2,2,3);
contour(x,y,real(f),nc);hold('on');plot([[x(1);0] [1;x(end)]],zeros(2),'k');
title('real');
subplot(2,2,4);
contour(x,y,imag(f),nc);hold('on');plot([[x(1);0] [1;x(end)]],zeros(2),'k');
title('imag');
