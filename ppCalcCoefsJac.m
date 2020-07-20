function [pCoefs,condNum] = ppCalcCoefsJac(a,psi0,nf,na,dzdx)

delta = 1/pi*acot(psi0);
gam = 1/2-1/pi*acot(psi0);

tcf = myJacobiNodes(nf,gam,1/2); % Generate forward collocation points
tca = myJacobiNodes(na,delta,gam); % Generate aft collocation points

t2fun = @(tVar) -1 + (tVar - 1)*(1+a)/(1-a);
t1fun = @(tVar)  1 + (tVar + 1)*(1-a)/(1+a);

xt1 = @(tVar) -1 + (tVar + 1)*(1+a)/2;
xt2 = @(tVar)  1 + (tVar - 1)*(1-a)/2;

%% Calculate the p-coefficients

% Set up forward section
jQ0f = weight(tcf,gam,-1/2).*myJacobiQ2(0,gam,-1/2,tcf);
jQAf = weight(tcf,gam,1/2).*myJacobiQ2(nf-2,gam,1/2,tcf);
jQBf = weight(t2fun(tcf),delta,gam).*myJacobiQ2(na-1,delta,gam,t2fun(tcf));

Ff = [-jQ0f,-jQAf,-jQBf];

% Set up aft section
jQ0a = weight(t1fun(tca),gam,-1/2).*myJacobiQ2(0,gam,-1/2,t1fun(tca));
jQAa = weight(t1fun(tca),gam, 1/2).*myJacobiQ2(nf-2,gam,1/2,t1fun(tca));
jPAa = weight(tca,delta,gam).*myJacobiP(na,na-1,delta,gam,tca);
jQBa = weight(tca,delta,gam).*myJacobiQ2(na-1,delta,gam,tca);

Fa = [-jQ0a, -jQAa, -jQBa + pi*psi0*jPAa];

% Now invert to get coefficients
F = [Ff;Fa];
condNum = cond(F);
rhsFun = 4*pi*[dzdx(xt1(tcf));dzdx(xt2(tca))];
pCoefs = F\rhsFun;

end