addpath('matlab2tikz/src')
%%
imageFolder = '../porous-jacobi/images/';


% Calculates the vorticity distribution along a partially porous aerofoil 
% using the Jacobi polynomial expansion method

% Set angle of attack and parabolic camber (only use alp=0 for now)
aoa = 1; bet = 0;

% Defines geometry and derivative
dzdx = @(xVar) -aoa - bet*xVar;

a = 0; 

psi0 = 0.1;

nT = 60; %total number of nodes;
nf = ceil(nT*(1+a)/2); na = floor(nT*(1-a)/2);

pCoefs = ppCalcCoefsJac(a,psi0,nf,na,dzdx);

% Extract forward and aft coefficients
pCoefsF = pCoefs(1:nf);
pCoefsA = pCoefs(nf+1:end);

% Set up plots
npf = 50;
npa = 50;

sGridF = .5*(1+sin(linspace(-pi,pi,npf).'/2));
sGridA = .5*(1+sin(linspace(-pi,pi,npa).'/2));
xPlotF = -1+(1+a)*sGridF;
xPlotA =  a+(1-a)*sGridA;

xPlot = [xPlotF;xPlotA];

pF = ppNumF(xPlotF,a,pCoefsF,psi0);
pA = ppNumA(xPlotA,a,pCoefsA,psi0);

numP = [pF;pA];

% Create plots    

figure(1)
plot(xPlot,abs(numP-ppExact(xPlot,a,psi0,aoa,bet)),'r');

figure(2)
plot(xPlot,abs(numP),'r');
hold on
plot(xPlot,abs(ppExact(xPlot,a,psi0,aoa,bet)),'bo')
hold off
axis([-1,1,0,10])


figure(3)
scatter(1:nf,abs(pCoefsF))
hold on
scatter(1:na,abs(pCoefsA))
hold off

close all
%% Calculate error
ne = 40;

np = 4;

%psis = linspace(1e-10,10,np);
%psis = 1;
%psis = 2*cot(flip(linspace(0,pi/2,np+2))); psis(1)=[]; psis(end) = [];
psis = cot(pi*flip([.1,.2,.3,.4]));
%psis = 0.1;
myErrorJac  = zeros(np,ne);
condNums= zeros(np,ne);

myErrorCheb1 = zeros(np,ne);
myErrorCheb2 = zeros(np,ne);

a=0;

nf = ceil((1:ne)*(1+a)/2)+2; na = floor((1:ne)*(1-a)/2)+1;
%nf = ceil((1:ne)/2+.5);
%na = 1:ne;
%nf = 50 + 0*(1:ne);
%na = 2:(ne+1);
%nf = 2:(ne+1);%ceil((1:ne)/2)+2; na = floor((1:ne)/2)+1;
%nf = (2:(ne+1)); + 0*(1:ne);
%na = (1:ne);

for l = 1:np

% Jacobi plots
    for k = 1:ne


        [pCoefsJac,condNum] = ppCalcCoefsJac(a,psis(l),nf(k),na(k),dzdx);
condNums(l,k) = condNum;



        % Extract forward and aft coefficients
        pCoefsFJac = pCoefsJac(1:nf(k));
        pCoefsAJac = pCoefsJac(nf(k)+1:end);

        integrandF = @(xVar) (ppExact(xVar,a,psis(l),aoa,bet)-ppNumF(xVar.',a,pCoefsFJac,psis(l)).').^2;
        integrandA = @(xVar) (ppExact(xVar,a,psis(l),aoa,bet)-ppNumA(xVar.',a,pCoefsAJac,psis(l)).').^2;
        xErrF = linspace(-1,a); xErrF(1) = []; xErrF(end) = [];
        xErrA = linspace(a,1); xErrA(1) = []; xErrA(end) = [];
        %maxErrorF = max(abs((ppExact(xErrF,a,psis(l),aoa,bet)-ppNumF(xErrF.',a,pCoefsFJac,psis(l)).')./ppNumF(xErrF.',a,pCoefsFJac,psis(l)).'));
        %maxErrorA = max(abs((ppExact(xErrA,a,psis(l),aoa,bet)-ppNumA(xErrA.',a,pCoefsAJac,psis(l)).')./ppNumA(xErrA.',a,pCoefsFJac,psis(l)).'));
        %maxErrorF = max(abs((ppExact(xErrF,a,psis(l),aoa,bet)-ppNumF(xErrF.',a,pCoefsFJac,psis(l)).')./ppExact(xErrF,a,psis(l),aoa,bet)));
%         %maxErrorA = max(abs((ppExact(xErrA,a,psis(l),aoa,bet)-ppNumA(xErrA.',a,pCoefsAJac,psis(l)).')./ppExact(xErrA,a,psis(l),aoa,bet)));
%         [maxErrorFA,locF] = max(abs((ppExact(xErrF,a,psis(l),aoa,bet)-ppNumF(xErrF.',a,pCoefsFJac,psis(l)).')));
%         [maxErrorAA,locA]=  max(abs((ppExact(xErrA,a,psis(l),aoa,bet)-ppNumA(xErrA.',a,pCoefsAJac,psis(l)).')));
%         maxErrorF = 0*maxErrorFA./abs(ppExact(xErrF(locF),a,psis(l),aoa,bet));
%         maxErrorA = maxErrorAA./abs(ppExact(xErrF(locA),a,psis(l),aoa,bet));
        
         myErrorF = abs(integral(integrandF,-1,a)./integral(@(xInt) ppExact(xInt,a,psis(l),aoa,bet).^2,-1,a));
         myErrorA = abs(integral(integrandA, a,1)./integral(@(xInt) ppExact(xInt,a,psis(l),aoa,bet).^2,a,1));
         myErrorJac(l,k) = myErrorF + myErrorA;
        %myErrorJac(l,k) = max([maxErrorF,maxErrorA]);
        %myErrorJac(l,k) = sqrt(integrandF(-.5));
        
            if k == 5
        figure(2)
        pF = ppNumF(xPlotF,a,pCoefsFJac,psis(l));
        pA = ppNumA(xPlotA,a,pCoefsAJac,psis(l));
        numP = [pF;pA];
        if l>1
            h = gca;
            set(h,'ColorOrderIndex',l);
        end
        scatter(xPlot(1:6:end),abs(numP(1:6:end)),'+','LineWidth',7);
        hold on
        h = gca;
        set(h,'ColorOrderIndex',l);
        plot(xPlot,abs(ppExact(xPlot,a,psis(l),aoa,bet)),'LineWidth',1.5)
        ylim([0,12])
        if l == np
                        xlabel('$x$','Interpreter','latex')
            if bet ==0
            ylabel('$\left|p / \alpha \right|$','Interpreter','Latex')
            cleanfigure;
                    %matlab2tikz([imageFolder,'pp-press-aoa','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

            else
            ylabel('$\left|p / \beta \right|$','Interpreter','Latex')
            cleanfigure;
                    %matlab2tikz([imageFolder,'pp-press-bet','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

            end
        end
    end
        
    end
% Now Cheb1 plots    
%     for k = 1:ne
% 
%         pCoefsCheb1 = ppCalcCoefsCheb(a,psis(l),nf(k),na(k),dzdx);
% 
%         % Extract forward and aft coefficients
%         pCoefsFCheb1 = pCoefsCheb1(1:nf(k));
%         pCoefsACheb1 = pCoefsCheb1(nf(k)+1:end);
% 
%         integrandF = @(xVar) (ppExact(xVar,a,psis(l),alph,bet)-ppNumFCheb1(xVar.',a,pCoefsFCheb1).').^2;
%         integrandA = @(xVar) (ppExact(xVar,a,psis(l),alph,bet)-ppNumACheb1(xVar.',a,pCoefsACheb1).').^2;
% 
%         myErrorF = abs(integral(integrandF,-1,a));
%         myErrorA = abs(integral(integrandA, a,1));
% 
%         myErrorCheb1(l,k) = myErrorF + myErrorA;
% 
%     end

   
% Now Cheb1 plots    
%     for k = 1:ne
% 
%         pCoefsCheb2 = ppCalcCoefsCheb2(a,psis(l),nf(k),na(k),dzdx);
% 
%         % Extract forward and aft coefficients
%         pCoefsFCheb2 = pCoefsCheb2(1:nf(k));
%         pCoefsACheb2 = pCoefsCheb2(nf(k)+1:end);
% 
%         integrandF = @(xVar) (ppExact(xVar,a,psis(l),alph,bet)-ppNumFCheb2(xVar.',a,pCoefsFCheb2).').^2;
%         integrandA = @(xVar) (ppExact(xVar,a,psis(l),alph,bet)-ppNumACheb2(xVar.',a,pCoefsACheb2).').^2;
% 
%         myErrorF = abs(integral(integrandF,-1,a));
%         myErrorA = abs(integral(integrandA, a,1));
% 
%         myErrorCheb2(l,k) = myErrorF + myErrorA;
% 
%     end

end
    hold off

%% Plot error

cols = lines;

figure(4)
for l = 1:np
    semilogy(abs(sqrt(myErrorJac(l,:))),'-','LineWidth',2,'Color',cols(l,:))
    hold on
 %    semilogy(sqrt(myErrorCheb1(l,:)),'o-','LineWidth',2,'Color',cols(l,:))
%     semilogy(sqrt(myErrorCheb2(l,:)),'x-','LineWidth',2,'Color',cols(l,:))
end

xlabel('$n$')
ylabel('$\epsilon_n$')

%ylim([1e-15,1])

hold off

cleanfigure;

if bet == 0;
%matlab2tikz([imageFolder,'pp-aoa-val','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
else
%matlab2tikz([imageFolder,'pp-cam-val','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
end

return

%% Plot chebyfun

figure(5)
plot(xPlotA,ppNumACheb(xPlotA,a,pCoefsA))
hold on
plot(xPlotA,ppExact(xPlotA,a,psi0,alph,bet))
plot(xPlotA,ppNumACheb(xPlotA,a,pCoefsAJac))

hold off

