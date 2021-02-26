% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';
% Set angle of attack and parabolic camber (only use alp=0 for now)
%beta0 = -1; beta1 = 1;

beta0 = 1; beta1 = 0;
z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) (beta1 + 0*xVar); struct.dzdx = dzdx;
%% Calculate the p-coefficients
clf

profile on
nk = 50;
%kVec = [eps,logspace(-2,3,nk)]; kVec(1) = [];
k1 = linspace(0,1.25,floor(2*nk/3)+1); k1(1) = [];
k2 = linspace(0,1,floor(nk/3)); k2(1) = [];
kVec = sort([log(1+k1).*k1,1./k2,logspace(1,2,10)]);
kPlot = [eps,.001,.01,.1,1,10];
kVec = sort([kVec,kPlot,100]); kVec = unique(kVec); %kVec(end) = []; % remove last entry for check
%kVec = 50;
nk = numel(kVec);
locs = [];
for l=1:numel(kPlot); locs = [locs,find(kVec==kPlot(l),1)];end

na = 4;
fullLift = zeros(nk,na,2);
ncLift = zeros(nk,na,2);
qsLift = zeros(nk,na,2);
fullCirc = zeros(nk,na,2);
ncCirc = zeros(nk,na,2);
qsCirc = zeros(nk,na,2);
%myncLift = zeros(nk,na);
%aVec = flip([-.5,-.25,0,.25,.5]);
psiFunCell = cell(na,1);
%ap = zeros(1,na);
%thet = zeros(1,na);
% Type loop: vary resistivity or effective density
for n = 1:2
    % Porosity loop
    for m = 1:na
        a = aVec(m);
% Frequency loop
        for j = 1:nk
            k = kVec(j);
            struct.N = round(10+20*sqrt(k)); % change 5 to 12.
            struct.k = k;

            if n == 1
                Psifun = @(x) 0.05*m*(1+x);
                rhoe = @(x) 1.5 + 0*x;
            elseif n ==2
                Psifun = @(x) 0.05*(1+x);
                rhoe = @(x) na+1-m+0*x;
            end

            Phifun = @(x) 1./Psifun(x);
            psiFunCell{m} = @(x) Psifun(x);
            struct.Phifun = Phifun;
            struct.rhoe = rhoe;

            solStruct = calculateUnsteadyCoefficients(struct); 

            [fLift,nLift,qLift] = lift(solStruct);
            [fCirc,nCirc,qCirc] = circulation(solStruct);
            qsLift(j,m,n)= qLift;
            ncLift(j,m,n) = nLift;
            fullLift(j,m,n) = fLift;
            qsCirc(j,m,n)= qCirc;
            fullCirc(j,m,n) = fCirc;
            disp([j,m,n]);
        end
    %struct.k = inf;
    %infSolStruct = calculateUnsteadyCoefficients(struct);
    %thet(m) = infSolStruct.Theta;
    end
end
profile off
%%
C = @(sigVar) besselk(1,1i*sigVar)./(besselk(0,1i*sigVar) + besselk(1,1i*sigVar));

%thetN = 1i*thet.*exp(1i*pi/4)*sqrt(2*pi);
theo = C(kVec);
numTheo = (fullLift - ncLift)./qsLift;

markerSize = 50;

for n = 1:2
    if n==1
    cols = flip(cmocean('matter',na+1));
    elseif n==2
    cols = flip(cmocean('speed',na+1));
    end

figure(1)
clf
plot(theo,'k','LineWidth',2)
hold on
scatter(real(theo(locs)),imag(theo(locs)),markerSize,'sk','MarkerFaceColor','k')
for m = flip(1:na)
plot(numTheo(:,m,n),'-','Color',cols(m,:),'LineWidth',1)
%plot(numTheo(:,m,n),'o','Color',cols(m,:),'LineWidth',1)
for j = 1:numel(kPlot)
scatter(real(numTheo(locs(j),m,n)),imag(numTheo(locs(j),m,n)),markerSize,'sk','MarkerFaceColor',cols(m,:))
end
end
%plot(theo,'k-','LineWidth',3)
grid on
hold on
hold off
%axis([-.1,1.1,-.35,.05])
axis([0,1.1,-.35,.05])
xlabel('$\Re[C(k)]$','Interpreter','latex')
ylabel('$\Im[C(k)]$','Interpreter','latex')

cleanfigure;
matlab2tikz([imageFolder,num2str(n),'theo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

% Absolute value plot
figure(2)
clf
for m = 1:na
semilogx(kVec,abs(numTheo(:,m,n)),'-','Color',cols(m,:),'LineWidth',1)
hold on
end
semilogx(kVec,abs(theo),'k','LineWidth',2);
grid on
hold off
ylabel('$|C(k)|$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
ylim([0,1])
xlim([1e-2,kVec(end)])

% ax2 = axes('Position',[.6 .6 .25 .25]);
% xPsi = linspace(-1,1);
% plot(xPsi+eps*1i,'k','LineWidth',2);
%     hold(ax2,'on')
% for m = flip(1:na)
%     plot(ax2,xPsi,.1*m*(1+xPsi), '-', 'LineWidth', 1,'Color',cols(m,:))
% end
% xlim(ax2,[-1,1])
% ylim(ax2,[-eps,Inf])
% xlabel(ax2,'$x$','Interpreter','latex')
% ylabel(ax2,'$1/\Phi$','Interpreter','latex')
% hold(ax2,'off')

cleanfigure;
matlab2tikz([imageFolder,num2str(n),'absTheo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

% Phase plot
figure(3)
clf
for m = 1:na
semilogx(kVec,angle(numTheo(:,m,n)),'-','Color',cols(m,:),'LineWidth',1)
hold on
end
semilogx(kVec,angle(theo),'k','LineWidth',2);
grid on
hold off
ylabel('$\angle C(k)$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
ylim([-.8,0])
xlim([1e-2,kVec(end)])

cleanfigure;
matlab2tikz([imageFolder,num2str(n),'phaseTheo.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
end