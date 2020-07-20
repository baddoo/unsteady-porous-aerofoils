% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi/images/';
% Set angle of attack and parabolic camber (only use alp=0 for now)
%beta0 = -1; beta1 = 1;

% z = @(xVar) 1+(xVar-.5);%+ xVar-.5+ 0*(xVar - .5);%beta0/2 + beta1*xVar; 
% struct.z = z;
%dzdx = @(xVar) +0*xVar.^2; struct.dzdx = dzdx;
beta0 = 1; beta1 = 0;

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) (beta1 + 0*xVar); struct.dzdx = dzdx;
%% Calculate the p-coefficients
clf

nN = 50;

na = 4;
errorMatJac = zeros(nN,na);
errorMatCheb = zeros(nN,na);
finalJacCoefs = zeros(nN+1,na);
finalChebCoefs = zeros(nN+1,na);

psiFunCell = {@(xVar) 1,...
             @(xVar) eps+(1+xVar)/2,...
             @(xVar) eps*.5*(1+xVar).^2,...
             @(xVar) .5+.25*sin(xVar)};
         
kVec = [.01,.1,1,5];
for m = 1:na
    struct.psiFun = psiFunCell{m};
    struct.k = kVec(m);
for N = 1:nN
    struct.N = N;

    jacStruct = struct; chebStruct = struct;
    jacStruct.type = 'full';
    jacSolStruct = calculateUnsteadyCoefficients(jacStruct); 
    chebStruct.type = 'cheb';
    chebSolStruct = calculateUnsteadyCoefficients(chebStruct); 
    if N>2
    errorMatJac(N,m) = sqrt(quadgk(@(xVar) (vortFun(xVar.',jacSolStruct).'-vortFun(xVar.',prevJacStruct).').^2,-1,1)./...
                  quadgk(@(xVar) vortFun(xVar.',jacSolStruct).'.^2,-1,1,'AbsTol',1e-12));
    errorMatCheb(N,m) = sqrt(quadgk(@(xVar) (vortFun(xVar.',chebSolStruct).'-vortFun(xVar.',prevChebStruct).').^2,-1,1)...
                  ./quadgk(@(xVar) vortFun(xVar.',chebSolStruct).'.^2,-1,1,'AbsTol',1e-12));
    end
    prevJacStruct = jacSolStruct;
    prevChebStruct = chebSolStruct;
    
    if N == nN
        finalJacCoefs(:,m) = jacSolStruct.coefs;
        finalChebCoefs(:,m)= chebSolStruct.coefs;
    end
end
end
return
%%
marker = ['s','o','x','p'];
cols = lines(4);
figure(1)
clf
for m = 1:na
    loglog(abs(errorMatCheb(:,m)),'s-','Color',cols(m,:))
    hold on
    loglog(abs(errorMatJac(:,m)),'x-','Color',cols(m,:))
end
xlim([3,nN])
hold off
figure(2)
clf
for m = 1:na
    semilogy(abs(finalChebCoefs(:,m)./finalChebCoefs(1,m)),'s-','Color',cols(m,:))
    hold on
    semilogy(abs(finalJacCoefs(:,m)./finalJacCoefs(1,m)),'x-','Color',cols(m,:))
end
xlim([1,nN])
hold off
return
%%
cols = (hot(floor(1.75*na)));
cols = flip(cols(1:na,:));
theo = C(kVec);
numTheo = (fullLift-ncLift)./qsLift;
theo(1) = 1;
numTheo(1,:) = 1;

% figure(1)
% loglog(kVec,abs(real(ncLift)),'-');
% hold on
% loglog(kVec,abs(imag(ncLift)),'--');
% loglog(kVec,abs(abs(ncLift)),'.-');
% hold off
% return
%numTheo = wakeLift./A;.*k

figure(1)
clf
for m = 1:na
semilogx(kVec,abs(numTheo(:,m)),'-','Color',cols(m,:),'LineWidth',2)
hold on
end
semilogx(kVec,abs(theo),'k','LineWidth',4);
grid on
hold off
ylabel('$|C(k)|$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
ylim([0,1.2])
xlim([1e-2,kVec(end)])

ax2 = axes('Position',[.6 .6 .25 .25]);
xPsi = linspace(-1,1);
plot(xPsi+eps*1i,'k','LineWidth',5);
    hold(ax2,'on')
for j = flip(1:na)
    plot(ax2,xPsi,psiFunCell{j}(xPsi), '-', 'LineWidth', 2,'Color',cols(j,:))
end
xlim([-1,1])
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$\psi$','Interpreter','latex')
hold(ax2,'off')

cleanfigure;
matlab2tikz([imageFolder,'absTheo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');


markerSize = 150;

figure(2)
clf
grid on
hold on
scatter(real(theo(locs)),imag(theo(locs)),markerSize,'sk','MarkerFaceColor','k')
for m = flip(1:na)
plot(numTheo(:,m),'-','Color',cols(m,:),'LineWidth',2)
for j = 1:numel(kPlot)
scatter(real(numTheo(locs(j),m)),imag(numTheo(locs(j),m)),markerSize,'sk','MarkerFaceColor',cols(m,:))
end
end
plot(theo,'k','LineWidth',4)
hold on
hold off
axis([-.1,1.1,-.35,.05])
xlabel('$\Re[C(k)]$','Interpreter','latex')
ylabel('$\Im[C(k)]$','Interpreter','latex')


%cleanfigure;
%matlab2tikz([imageFolder,'theo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
