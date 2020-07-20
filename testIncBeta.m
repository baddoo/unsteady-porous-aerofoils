z = linspace(0,1); z(1) = []; z(end) = [];

a = .75; b = -.5;

tic
%f1 = betainc(z,b,a)*beta(b,a);
f1 = integral(@(t) heaviside(z-t).*t.^(b-1).*(1-t).^(a-1),0,1,'ArrayValued',true)
toc
tic
f2 = inbeta(z,b,a,100);
toc

err = norm(f1-f2,'inf')

plot(f1,'LineWidth',5)
hold on
plot(f2,'--','LineWidth',5)
hold off