

teZero = .25;
leSing = .75;
gamCoefs = [2,1,1,3].';
k=1;
vortFunInt = @(zVar) vortFun(zVar.',gamCoefs, teZero,leSing,k).';
%presFun = @(zVar) pFun2(zVar.',gamCoefs, teZero,leSing,k).';

numInt = @(xVar) vortFun(xVar.',gamCoefs, teZero,leSing,k).'  + 1i*k*integral(vortFunInt,-1,xVar,'ArrayValued',true);

p = @(xVar) presFun(xVar.',gamCoefs,teZero,leSing,k).';
%% Plots
xPlot = sin(pi/2*linspace(-1,1)); xPlot(1) = []; xPlot(end) = [];

plotInt = arrayfun(numInt,xPlot);
figure(1) 
clf
plot(xPlot,real(plotInt-p(xPlot)));
hold on
%plot(xPlot,real(p(xPlot)),'.');
hold off