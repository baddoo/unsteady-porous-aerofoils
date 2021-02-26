% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';

%% Calculate the p-coefficients
clf

nPsi = 1;
psi = linspace(0,2,nPsi);
N = 30;
struct.N = N; 
struct.type = 'full';
nj = 5;
PSIVec = linspace(0,5,nj);
xp = cos(flip(linspace(0,pi,1e2))'); xp(1) = []; xp(end) = [];
xp = sort([xp;1-logspace(-5,-1,10)']);
pMat = zeros(nj,numel(xp));
psiFunCell = cell(1,nj);

profile on
for m = 1:2
    
% Set angle of attack and parabolic camber (only use alp=0 for now)
if m ==1; beta0 = 1; beta1 = beta0/2; elseif m==2; beta0 = 1; beta1 = 0; end
% m ==1 is pure pitching, m==2 is pure heaving.
z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) beta1 + 0*xVar; struct.dzdx = dzdx;

for j = 1:nj

if m ==1
    struct.k = .5;
elseif m==2
    struct.k = .1;
end

if m == 1
psiFun = @(x) eps+.1*PSIVec(j).*(x + 1);
elseif m==2
psiFun = @(x) eps+0.05*PSIVec(j).*(x + 1).^2;
end

rhoe = @(x) 1.2 + 0*x;
struct.Phifun = @(x) 1./psiFun(x);
struct.rhoe = rhoe;

psiFunCell{j} = psiFun;
tic
solStruct = calculateUnsteadyCoefficients(struct); 
toc
pMat(j,:) = presFun(xp,solStruct)';

end

profile off
pMat(:,end) = 0;

%% Plots
cols = hot(ceil(1.5*nj));

figure(m)
clf
hold on
plot(xp,abs(pMat(1,:)),'k-','LineWidth',2);
for j = 2:nj
plot(xp,abs(pMat(j,:)),'-','LineWidth',1,'Color',cols(j,:));
ax1 = gca;
end
hold off
xlim([-1,1]);
grid on

if m==1; ylim([0,5]); elseif m ==2; ylim([0,.5]); end

xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$|\Delta p|$','Interpreter','latex')


ax2 = axes('Position',[.6 .6 .25 .25]);
     hold(ax2,'on');
plot(ax2,xp,real(psiFunCell{1}(xp)), 'k-', 'LineWidth', 2)
for j = 2:nj
    plot(ax2,xp,real(psiFunCell{j}(xp)), '-', 'LineWidth', 1,'Color',cols(j,:))
    %semilogy(ax2,xFMat(j,:),psiFunCell{j}(xFMat(j,:)), '-', 'LineWidth', 2,'Color',cols(j,:))
end
xlim([-1,1])
grid minor
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$1/\Phi$','Interpreter','latex')
hold off
cleanfigure;
matlab2tikz([imageFolder,'largePorCont',num2str(m),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
                 'extraaxisoptions','axis on top=false');

end

