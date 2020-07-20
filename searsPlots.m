% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
%imageFolder = '../unsteady-jacobi/images/';
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';
% Set angle of attack and parabolic camber (only use alp=0 for now)

z = @(xVar) 0*xVar; struct.z = z;

%% Calculate the p-coefficients
clf

profile on
tic

nk = 80;
kPlot = [eps,.1,1,2,5,10];
kVec = [logspace(-2,1,nk/2),linspace(1,10,nk/2)];
kVec = sort([kVec,kPlot]); kVec = unique(kVec);

nP = 4;
psiVec = [.1,.2,.5,1.5];%logspace(-1,1,nP);
LqsPP= zeros(1,nP);
newnk = numel(kVec);
fullLift = zeros(newnk,nP);
qsLift = zeros(newnk,nP);

ppLift = zeros(newnk,nP);
psiFunCellPP = cell(1,nP);
psiFunCell = cell(1,nP);
a = -.5; b = .5;
aVec = flip([-.6,-.3,0,.3]);
bVec = aVec+.3;
psi0 = 1;
na = 15; nf = 20;
for m = 1:nP
    junction = aVec(m)+.15;
    Psifun = @(xVar) eps + .25*psi0*(heaviside((xVar-aVec(m)).*(bVec(m)-xVar)).*(xVar-aVec(m))/(bVec(m)-aVec(m))...
                                     +heaviside(xVar-bVec(m)));
    Phifun = @(x) 1./Psifun(x);
    rhoe = @(x) 1 + 0*x;
    psiFunCellPP{m} = @(xVar) eps + psi0*heaviside(xVar-junction);
    psiFunCell{m} = @(x) Psifun(x);
    struct.Phifun = Phifun;
    struct.rhoe = rhoe;

    %solStructQS = calculateUnsteadyCoefficients(qsStruct); 
    
%     qsStructPP = qsStruct;
%     qsStructPP.na = na; qsStructPP.nf = nf;
%     qsStructPP.psiFun = psiFunCellPP{m};
%     qsStructPP.junction = junction;
%     solStructQSPP = calculateUnsteadyCoefficientsDiscont(qsStructPP); 
      %qsLift(m) = integral(@(xVar) presFun(xVar',solStructQS).',-1,1);
%     LqsPP(m) = quadgk(@(xInt) presFUnsteady(xInt.',solStructQSPP).',-1,junction) ...
%             +quadgk(@(xInt) presAUnsteady(xInt.',solStructQSPP).', junction,1);
dzdx = @(xVar) -1; struct.dzdx = dzdx;
structQS = struct;
structQS.k = eps;

    for j = 1:newnk 
        k = kVec(j);
        struct.k = k;
        struct.N = floor(10+12*sqrt(k));
        dzdx = @(xVar) -exp(-1i*k*xVar); struct.dzdx = dzdx;
         structPP=struct; 
         structPP.na = na; structPP.nf = nf;
        structPP.psiFun = @(x) 1./psiFunCellPP{m}(x);
        structPP.junction = junction;
        solStruct = calculateUnsteadyCoefficients(struct);   
        [fLift,~,qLift] = lift(solStruct);
        fullLift(j,m) = fLift;
        qsLift(j,m) = qLift;
        
%         solStructPP = calculateUnsteadyCoefficientsDiscont(structPP); 
%         ppLift(j,m) = quadgk(@(xInt) presFUnsteady(xInt.',solStructPP).',-1,junction) ...
%                      + quadgk(@(xInt) presAUnsteady(xInt.',solStructPP).', junction,1);
disp([m,j])
    end
end
toc
profile off

S = @(kVar) -1i./(kVar.*(besselk(0,1i*kVar) + besselk(1,1i*kVar)));

%%
sears = S(kVec);
numSears = fullLift/(-4*pi);%/Lqs;
numSearsNorm = fullLift./qsLift;
%numSearsPP = ppLift/(-4*pi);
%numSearsNormPP = ppLift./LqsPP;

cols = hot(floor(1*(nP+2)));
markerSize = 150;
markers = 'sssssssss';

locs = [];
for l=1:numel(kPlot); locs = [locs,find(kVec==kPlot(l),1)];end

figure(1)
clf
plot(sears,'k','LineWidth',4);
ax1 = gca;
grid on
hold on
for j = 1:numel(kPlot)
scatter(real(sears(locs(j))),imag(sears(locs(j))),markerSize,[markers(j),'k'],'MarkerFaceColor','k')
end
for m = 1:nP
plot(numSearsNorm(:,m),'LineWidth',2,'Color',cols(m,:));
%plot(numSearsNormPP(:,m),'--','LineWidth',2,'Color',cols(m,:));
for j = 1:numel(kPlot)
scatter(real(numSearsNorm(locs(j),m)),imag(numSearsNorm(locs(j),m)),markerSize,[markers(j),'k'],...
    'MarkerFaceColor',cols(m,:))
%scatter(real(numSearsNormPP(locs(j),m)),imag(numSearsNormPP(locs(j),m)),markerSize,[markers(j),'k'],...
%    'MarkerFaceColor',cols(m,:))
end
end
hold off
axis([-.3,1.1,-.2,.3])
xlabel(ax1,'$\Re[\mathcal{S}]$','Interpreter','latex')
ylabel(ax1,'$\Im[\mathcal{S}]$','Interpreter','latex')

%cleanfigure;
%matlab2tikz([imageFolder,'Sears.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');


figure(2)
clf
ax1 = gca;
plot(sears,'k','LineWidth',4);
grid on
hold on
for j = 1:numel(kPlot)
scatter(real(sears(locs(j))),imag(sears(locs(j))),markerSize,[markers(j),'k'],'MarkerFaceColor','k')
end
for m = flip(1:nP)
plot(numSears(:,m),'-','LineWidth',2,'Color',cols(m,:));
%plot(numSearsPP(:,m),'--','LineWidth',2,'Color',cols(m,:));
for j = 1:numel(kPlot)
scatter(real(numSears(locs(j),m)),imag(numSears(locs(j),m)),markerSize,[markers(j),'k'],'MarkerFaceColor',cols(m,:))
%scatter(real(numSearsPP(locs(j),m)),imag(numSearsPP(locs(j),m)),markerSize,[markers(j),'k'],'MarkerFaceColor',cols(m,:))

end
end
hold off
xlabel(ax1,'$\Re[\mathcal{L}]$','Interpreter','latex')
ylabel(ax1,'$\Im[\mathcal{L}]$','Interpreter','latex')
axis([-.3,1.1,-.2,.3])

ax2 = axes('Position',[.6 .6 .25 .25]);
xPsi = linspace(-1,1);
    hold(ax2,'on')
for j = flip(1:nP)
    plot(ax2,xPsi,psiFunCell{j}(xPsi), '-', 'LineWidth', 2,'Color',cols(j,:))
    %plot(ax2,xPsi,psiFunCellPP{j}(xPsi), '--', 'LineWidth', 2,'Color',cols(j,:))
end
plot(xPsi+eps*1i,'k','LineWidth',4);

xlim([-1,1])
ylim([0,.31])
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$1/\Phi$','Interpreter','latex')
hold(ax2,'off')

%cleanfigure;
%matlab2tikz([imageFolder,'LiftVar.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
