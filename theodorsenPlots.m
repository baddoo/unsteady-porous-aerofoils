% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';
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

profile on
nk = 20;
%kVec = [eps,logspace(-2,3,nk)]; kVec(1) = [];
kVec = logspace(-2.5,2.5,nk);
kPlot = [eps,.001,.01,.1,1,10];
kVec = sort([kVec,kPlot]); kVec = unique(kVec); %kVec(end) = []; % remove last entry for check
%kVec = 50;
nk = numel(kVec);
locs = [];
for l=1:numel(kPlot); locs = [locs,find(kVec==kPlot(l),1)];end

na = 4;
fullLift = zeros(nk,na);
ncLift = zeros(nk,na);
qsLift = zeros(nk,na);
fullCirc = zeros(nk,na);
ncCirc = zeros(nk,na);
qsCirc = zeros(nk,na);
%myncLift = zeros(nk,na);
aVec = flip([-.5,-.25,0,.25,.5]);
psiFunCell = cell(na,1);
ap = zeros(1,na);
thet = zeros(1,na);
for m = 1:na
a = aVec(m);
for j = 1:nk
    k = kVec(j);
    struct.N = round(10+20*sqrt(k)); % change 5 to 12.
    struct.k = k;
    
    Psifun = @(x) 0.1*m*(1+x);
    Phifun = @(x) 1./Psifun(x);
    psiFunCell{m} = @(x) Psifun(x);
    rhoe = @(x) 1 + 0*x;
    struct.Phifun = Phifun;
    struct.rhoe = rhoe;
    
    solStruct = calculateUnsteadyCoefficients(struct); 

    [fLift,nLift,qLift] = lift(solStruct);
    [fCirc,nCirc,qCirc] = circulation(solStruct);
    qsLift(j,m)= qLift;
    ncLift(j,m) = nLift;
    fullLift(j,m) = fLift;
    qsCirc(j,m)= qCirc;
    %ncCirc(j,m) = nCirc;
    fullCirc(j,m) = fCirc;
    %wakeLift(j,m) = wLift;
    disp([j,m])
end
struct.k = inf;
infSolStruct = calculateUnsteadyCoefficients(struct);
thet(m) = infSolStruct.Theta;
end

C = @(sigVar) besselk(1,1i*sigVar)./(besselk(0,1i*sigVar) + besselk(1,1i*sigVar));

profile off
%%
thetN = 1i*thet.*exp(1i*pi/4)*sqrt(2*pi);
cols = (hot(floor(1.5*(na+1))));
theo = C(kVec);
numTheo = (fullLift - ncLift)./qsLift;

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
for m = flip(1:na)
    plot(ax2,xPsi,.1*m*(1+xPsi), '-', 'LineWidth', 2,'Color',cols(m,:))
end
xlim(ax2,[-1,1])
ylim(ax2,[0,.6])
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$1/\Phi$','Interpreter','latex')
hold(ax2,'off')

%cleanfigure;
%matlab2tikz([imageFolder,'absTheo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

markerSize = 150;

figure(2)
clf
plot(theo,'k','LineWidth',7)
hold on
scatter(real(theo(locs)),imag(theo(locs)),markerSize,'sk','MarkerFaceColor','k')
for m = flip(1:na)
plot(numTheo(:,m),'-','Color',cols(m,:),'LineWidth',5)
for j = 1:numel(kPlot)
scatter(real(numTheo(locs(j),m)),imag(numTheo(locs(j),m)),markerSize,'sk','MarkerFaceColor',cols(m,:))
end
end
%plot(theo,'k-','LineWidth',3)
grid on
hold on
hold off
%axis([-.1,1.1,-.35,.05])
axis([0,1.2,-.4,.1])
xlabel('$\Re[C(k)]$','Interpreter','latex')
ylabel('$\Im[C(k)]$','Interpreter','latex')

%cleanfigure;
%matlab2tikz([imageFolder,'theo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
