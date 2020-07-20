addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';

INT = 'interpreter';
LW = 'LineWidth';
CL = 'Color';

na = 5;
PHI = logspace(-2,2,1e3);
rhoe = 1;%+linspace(0,.5,na)';
RHO = linspace(1,100,1e3);
k1 =(0:na-1)';
k2 = .25*(0:na-1)';

psi1 = -4./(2i*k1*rhoe + PHI);
psi2 = -4./(2i*k2*RHO + 1);

cols = cool(na);

arccot = @(x) 1i/2*(log((1 - 1i./x)) - log(1*(1 + 1i./x)));
alph1 = 1/pi*arccot(psi1);
alph2 = 1/pi*arccot(psi2);


figure(1)
for j = 1:na
semilogx(PHI,abs(real(alph1(j,:))),CL,cols(j,:),LW, 2)
hold on
semilogx(PHI,abs(imag(alph1(j,:))),'--',CL,cols(j,:), LW, 2)
end
hold off

set(gca,'XDir','reverse');
ylabel('$\frac{1}{\pi} \cot^{-1}(\psi)$',INT,'Latex')
xlabel('$\Phi$',INT,'Latex')
ylim([0,.5])
grid on
grid off
grid on

cleanfigure;
matlab2tikz([imageFolder,'alphVarPsi.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,...
            'extratikzpictureoptions','trim axis left, trim axis right');


figure(2)
for j = 1:na
semilogx(RHO,abs(real(alph2(j,:))),CL,cols(j,:),LW, 2)
hold on
semilogx(RHO,abs(imag(alph2(j,:))),'--',CL,cols(j,:), LW, 2)
end
hold off
ylabel('$\frac{1}{\pi} \cot^{-1}(\psi)$',INT,'Latex')
xlabel('$\rho_e$',INT,'Latex')
ylim([0,.5])

grid on
grid off
grid on
cleanfigure;
matlab2tikz([imageFolder,'alphVarRho.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,...
            'extratikzpictureoptions','trim axis left, trim axis right');

