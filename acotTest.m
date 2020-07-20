x = linspace(-2,2)
z = x + 1i*x';

a1 = @(x) acot(x);
a2 = @(x) 1i/2*(log(-(x-1i)./(x+1i))-1i*pi);

figure(1)
pcolor(real(z),imag(z),real(a1(z)))
shading interp
colorbar
caxis([-pi,pi])

figure(2)
pcolor(real(z),imag(z),real(a2(z)))
shading interp
colorbar
caxis([-pi,pi])
