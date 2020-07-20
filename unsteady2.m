% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion

% Set angle of attack and parabolic camber (only use alp=0 for now)
beta0 = .5; beta1 = .5;

freq = 1;
k=2*pi*freq;
U = 1/freq;
% Defines geometry and derivative
struct.z = @(xVar) beta0/2 + beta1*xVar;
struct.dzdx = @(xVar) beta1;
%struct.wa = @(xVar) 1i*k*z(xVar)+dzdx(xVar);
struct.k = k;
%% Calculate the p-coefficients
nVar = 10;
clf

psiV = logspace(-3,2,nVar);
psiFunCell = cell(1,nVar);
for l = 1:nVar
    %psiFunCell{l} = @(xVar) psiV(l) + 0*xVar;
    psiFunCell{l} = @(xVar) psiV(l)*(1+xVar)/2;
    psiFunCell{l} = @(xVar)  psiV(l).*(1+xVar).^2/4;
    %psiFunCell{l} = @(xVar)  psiV(l).*sin(5*(1+xVar));
end

fullCell = cell(1,nVar);

struct.N = 100;

profile on
tic
for l = 1:nVar
   
struct.psiFun = psiFunCell{l};
solStruct = calculateUnsteadyCoefficients(struct);   
fullCell{l} = @(zVar) presFun(zVar.',solStruct).';

end

toc
profile off

%% Plots
nP = 500;
xPlot = sin(pi/2*linspace(-1,1,nP+2)); xPlot(1) = []; xPlot(end) = [];
cols = (hot(ceil(1.5*nVar)));

figure(1)
clf; cla

ax1 = gca;
plot(ax1,xPlot, abs(fullCell{1}(xPlot)),'-' ,'LineWidth',4,'Color','k')    
hold on
for l0 = 2:nVar
    %ncvals = ncFunCell{l0}(xPlot);
    %cvals = cFunCell{l0}(xPlot);
    %semilogy(ax1,xPlot, abs(ncvals + wakeStrength(l0)*cvals),'-' ,'LineWidth',4,'Color',cols(l0,:))  
    
    plot(ax1,xPlot, abs(fullCell{l0}(xPlot)),'-' ,'LineWidth',1,'Color',cols(l0,:))    

    hold(ax1,'on');
end

set(ax1,'ylim',[0,100])
set(ax1,'xlim',[-1,1]);

% %xlabel(ax1,'$x$','Interpreter','latex')
% %ylabel(ax1,'$\gamma$','Interpreter','latex')
% 
% ax2 = axes('Position',[.2 .2 .25 .25]);
% for l0 = 1:nVar
%     plot(ax2,xPlot,psiFunCell{l0}(xPlot) , '-', 'LineWidth', 2,'Color',cols(l0,:))
%      hold(ax2,'on');
% end
% xlabel(ax2,'$x$','Interpreter','latex')
% ylabel(ax2,'$\psi$','Interpreter','latex')
% hold off

sig = 2*pi/U;
C = @(sigVar) besselk(1,1i*sigVar)/(besselk(0,1i*sigVar) + besselk(1,1i*sigVar));

a = zeros(1,3);
a(1) = -2i*pi*U*C(sig)*beta0 + 2i*pi*U*(1-C(sig))*beta1 - 2*C(sig)*U^2*beta1;
a(2) = 2*pi^2*beta0 - 4i*pi*U*beta1;
a(3) = pi^2*beta1;

pMoore = @(xVar) a(1)*sqrt((1-xVar)./(1+xVar)) ...
          + a(2)*2*sqrt(1-xVar.^2) ...
          + a(3)*4*xVar.*sqrt(1-xVar.^2);
plot(xPlot,abs(pMoore(xPlot)/U.^2),'g')
hold off
%%
return
figure(2) 
error = zeros(nVar-1,1);
for l0 = 1:(nVar-1)
     intFun = @(xVar) abs(fullCell{l0+1}(xVar) - fullCell{l0}(xVar)).^2;
     norFun = @(xVar) abs(fullCell{l0+1}(xVar)).^2;
     error(l0)= sqrt(quadgk(intFun,-1,1)./quadgk(norFun,-1,1));
    %error(l0) = sqrt(max(intFun(xPlot))./max(norFun(xPlot)));
end

length = 1:nVar-1;
loglog(length,abs(error));
