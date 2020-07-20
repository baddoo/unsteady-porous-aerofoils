
k = 50;2*pi;

beta0 = 1; beta1 = 0;

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) (beta1 + 0*xVar); struct.dzdx = dzdx;
struct.k = k;
struct.N = 300;
struct.psiFun = @(xVar) eps+ 0*(1+xVar)/2;
solStruct = calculateUnsteadyCoefficients(struct);   

sig = k;
C = @(sigVar) besselk(1,1i*sigVar)/(besselk(0,1i*sigVar) + besselk(1,1i*sigVar));

U = 2*pi/k;

a = zeros(1,3);
a(1) = -2i*pi*C(sig)*U*beta0 + 2i*pi*U*(1-C(sig))*beta1 - 2*C(sig)*U^2*beta1;
a(2) = 2*pi^2*beta0 - 4i*pi*U*beta1;
a(3) = pi^2*beta1;

pMoore = @(xVar) a(1)*sqrt((1-xVar)./(1+xVar)) ...
          + a(2)*2*sqrt(1-xVar.^2) ...
          + a(3)*4*xVar.*sqrt(1-xVar.^2);
      
%const = 2*pi^2;      
const = -2*pi^2/k^2;

lift = a(1)*2^(1-.5+.5)*beta(1-.5,1+.5) + 2*a(2)*2^(1+.5+.5)*beta(1+.5,1+.5)
liftNum = integral(@(xInt) pMoore(xInt),-1,1)
liftNum = integral(@(xInt) presFun(xInt',solStruct).',-1,1)*const
lift - liftNum
%% Plots
figure(1)
xPlot = sin(linspace(-1,1)*pi/2); xPlot(1)=[]; xPlot(end) = [];
moorePres = pMoore(xPlot);
myPres = presFun(xPlot',solStruct).'*const;
semilogy(xPlot,abs((moorePres-myPres)./moorePres),'LineWidth',4)
hold on
%semilogy(xPlot,imag()*const,'--','LineWidth',2)
hold off
%ylim([.01,100])