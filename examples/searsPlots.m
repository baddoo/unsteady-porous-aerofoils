% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
%imageFolder = '../unsteady-jacobi/images/';
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';
z = @(xVar) 0*xVar; struct.z = z;

%% Calculate the p-coefficients
clf

profile on
nk = 160;
kPlot = [eps,.1,1,2,5,10];
kVec = [logspace(-2,2,180),linspace(1,10,nk)];
kVec = sort([kVec,kPlot]); kVec = unique(kVec);

nP = 4;
LqsPP= zeros(1,nP);
newnk = numel(kVec);
fullLift = zeros(newnk,nP,2);
qsLift = zeros(newnk,nP,2);

ppLift = zeros(newnk,nP);
psiFunCellPP = cell(1,nP);
psiFunCell = cell(1,nP);
a = -.5; b = .5;
aVec = flip([-.5,-.25,0,.25]);bVec = aVec+.8;

psi0 = .5;
tic
na = 15; nf = 20;
for n = 1:2
for m = 1:nP
    junction = aVec(m)+.15;
    if n==1
        Psifun = @(xVar) eps + psi0*(heaviside((xVar-aVec(m)).*(bVec(m)-xVar)).*(xVar-aVec(m))/(bVec(m)-aVec(m))...
                                     +heaviside(xVar-bVec(m)));
        rhoe = @(x) 1.5+ 0*x;
    elseif n==2
        sifun = @(xVar) eps + 0.5*(heaviside((xVar-aVec(1)).*(bVec(1)-xVar)).*(xVar-aVec(1))/(bVec(1)-aVec(1))...
                                     +heaviside(xVar-bVec(1)));
        rhoe = @(x) nP+1-m + 0*x;
    end
    Phifun = @(x) 1./Psifun(x);
    psiFunCell{m} = @(x) Psifun(x);
    struct.Phifun = Phifun;
    struct.rhoe = rhoe;
dzdx = @(xVar) -1; struct.dzdx = dzdx;
structQS = struct;
structQS.k = eps;

    for j = 1:newnk 
        k = kVec(j);
        struct.k = k;
        struct.N = floor(25+20*sqrt(k));
        dzdx = @(xVar) -exp(-1i*k*xVar); struct.dzdx = dzdx;
        solStruct = calculateUnsteadyCoefficients(struct);   
        [fLift,~,qLift] = lift(solStruct);
        fullLift(j,m,n) = fLift;
        qsLift(j,m,n) = qLift;
        disp([m,j,n])
    end
end
end
   toc

profile off

S = @(kVar) -1i./(kVar.*(besselk(0,1i*kVar) + besselk(1,1i*kVar)));

%%
sears = S(kVec);
numSears = fullLift/(-4*pi);%/Lqs;
%numSearsNorm = fullLift./qsLift;
numSearsNorm = fullLift./fullLift(1,:,:);
%numSearsPP = ppLift/(-4*pi);
%numSearsNormPP = ppLift./LqsPP;

markerSize = 50;
markers = 'sssssssss';

locs = [];
for l=1:numel(kPlot); locs = [locs,find(kVec==kPlot(l),1)];end

locsP = 1:find(kVec==10);

for n = 1:2
    
    if n==1
    cols = flip(cmocean('matter',nP+1));
    elseif n==2
    cols = flip(cmocean('speed',nP+1));
    end    
% Sears
figure(1)
clf
plot(sears(locsP),'k','LineWidth',2);
ax1 = gca;
grid on
hold on
for j = 1:numel(kPlot)
scatter(real(sears(locs(j))),imag(sears(locs(j))),markerSize,[markers(j),'k'],'MarkerFaceColor','k')
end
for m = 1:nP
plot(numSearsNorm(locsP,m,n),'LineWidth',1,'Color',cols(m,:));
for j = 1:numel(kPlot)
scatter(real(numSearsNorm(locs(j),m,n)),imag(numSearsNorm(locs(j),m,n)),markerSize,[markers(j),'k'],...
    'MarkerFaceColor',cols(m,:))
end
end
hold off
axis([-.35,1.1,-.3,.45])
xlabel(ax1,'$\Re[S(k)]$','Interpreter','latex')
ylabel(ax1,'$\Im[S(k)]$','Interpreter','latex')

cleanfigure;
matlab2tikz([imageFolder,num2str(n),'Sears.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
% Absolute value
% Sears
figure(2)
clf
semilogx(kVec,abs(sears),'k','LineWidth',2);
ax1 = gca;
grid on
hold on
for m = 1:nP
semilogx(kVec,abs(numSearsNorm(:,m,n)),'LineWidth',1,'Color',cols(m,:));
end
hold off
%axis([-.3,1.1,-.2,.3])
ylabel(ax1,'$|S(k)|$','Interpreter','latex')
xlabel(ax1,'$k$','Interpreter','latex')
xlim([1e-2,1e2])
cleanfigure;
matlab2tikz([imageFolder,num2str(n),'absSears.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

figure(3)
% Phase
clf
semilogx(kVec,angle(exp(-1i*kVec).*sears),'k','LineWidth',2);
ax1 = gca;
grid on
hold on
for m = 1:nP
semilogx(kVec,angle(exp(-1i*kVec').*numSearsNorm(:,m,n)),'LineWidth',1,'Color',cols(m,:));
end
hold off
%axis([-.3,1.1,-.2,.3])
ylabel(ax1,'$\angle \e^{-\textrm{i} k}S(k)$','Interpreter','latex')
xlabel(ax1,'$k$','Interpreter','latex')
xlim([1e-2,1e2])

cleanfigure;
matlab2tikz([imageFolder,num2str(n),'phaseSears.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
end

return
% figure(2)
% clf
% ax1 = gca;
% plot(sears,'k','LineWidth',4);
% grid on
% hold on
% for j = 1:numel(kPlot)
% scatter(real(sears(locs(j))),imag(sears(locs(j))),markerSize,[markers(j),'k'],'MarkerFaceColor','k')
% end
% for m = flip(1:nP)
% plot(numSears(:,m),'-','LineWidth',2,'Color',cols(m,:));
% %plot(numSearsPP(:,m),'--','LineWidth',2,'Color',cols(m,:));
% for j = 1:numel(kPlot)
% scatter(real(numSears(locs(j),m)),imag(numSears(locs(j),m)),markerSize,[markers(j),'k'],'MarkerFaceColor',cols(m,:))
% %scatter(real(numSearsPP(locs(j),m)),imag(numSearsPP(locs(j),m)),markerSize,[markers(j),'k'],'MarkerFaceColor',cols(m,:))
% 
% end
%end
% hold off
% xlabel(ax1,'$\Re[\mathcal{L}]$','Interpreter','latex')
% ylabel(ax1,'$\Im[\mathcal{L}]$','Interpreter','latex')
% axis([-.3,1.1,-.2,.3])
% 
% ax2 = axes('Position',[.6 .6 .25 .25]);
xPsi = linspace(-1,1);
    %hold(ax2,'on')
for j = flip(1:nP)
            psiFun = @(xVar) eps + psi0*(heaviside((xVar-aVec(j)).*(bVec(j)-xVar)).*(xVar-aVec(j))/(bVec(j)-aVec(j))...
                                     +heaviside(xVar-bVec(j)));
    plot(xPsi,psiFun(xPsi), '-', 'LineWidth', 1,'Color',cols(j,:))
    hold on
    %plot(ax2,xPsi,psiFunCellPP{j}(xPsi), '--', 'LineWidth', 2,'Color',cols(j,:))
end
plot(xPsi+eps*1i,'k','LineWidth',4);
grid on
xlim([-1,1])
ylim([0,0.6])
xlabel('$x$','Interpreter','latex')
ylabel('$1/\Phi$','Interpreter','latex')
hold off
cleanfigure;
matlab2tikz([imageFolder,'searsPorosity.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
