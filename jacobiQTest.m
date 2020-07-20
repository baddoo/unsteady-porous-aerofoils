nx = 1e1; xInt = linspace(-1,1,nx+2).'; xInt(1)=[];xInt(end)=[];

alph = 0.35; bet = .35;
N = 20;
trans = 32;
Q = weight(xInt+trans,alph,bet).*myJacobiQ2(N-1,alph,bet,xInt+trans);
myN = 15;

%%

integrandU = @(xVar) weight(xVar,alph,bet).*jacobiP(myN,alph,bet,xVar);
%regIntegrand = @(xVar) (integrandU(xVar)-integrandU(xInt))./(xVar-xInt);
regIntegrand = @(xVar) integrandU(xVar)./(xVar-(xInt+trans));
regInt = integral(regIntegrand,-1,1,'ArrayValued',true);
%intFun = regInt + integrandU(xInt).*log((1-xInt)./(1+xInt));
intFun = regInt;

%%
plot(xInt,intFun,'k','LineWidth',2)
hold on
plot(xInt,Q(:,myN+1),'ro')
hold off