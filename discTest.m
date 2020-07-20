% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi/images/';

% Set angle of attack and parabolic camber (only use alp=0 for now)
beta0 = 1; beta1 = 0;

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) beta1 + 0*xVar; struct.dzdx = dzdx;

%% Calculate the p-coefficients
clf

nVar = 100;

solCell = cell(1,nVar);

nPsi = 1;
psi = linspace(0,2,nPsi);
struct.k = 1;
struct.nf =2; struct.na = 2;
nj = 5;
psiVec = logspace(-1,2,nj);
nJunct = 5;
aVec = linspace(-1,1,nJunct+2);aVec(1) = []; aVec(end) = [];
nx = 100;
xAMat = zeros(nj,nx);xFMat = zeros(nj,nx);
pAMat = zeros(nj,nx);pFMat = zeros(nj,nx);
vAMat = zeros(nj,nx);vFMat = zeros(nj,nx);
liftMat = zeros(nJunct,1);
psiFunCell = cell(1,nJunct);

for m = 1:2
for j = 1:nJunct

a = aVec(j);
struct.junction = a;
if m ==1
struct.psiFun = @(xVar) .1*heaviside(xVar-a)./(1-a);
elseif m==2
struct.psiFun = @(xVar) heaviside(xVar-a).*(1/(1-a)+.5*(xVar-a)/(1-a).^2);
end
psiFunCell{j} = struct.psiFun;
solStruct = calculateUnsteadyCoefficientsDiscont(struct);   

xf = -1+(1+a)*(1+sin(linspace(-1,1,nx)*pi/2)).'/2;
xa = a+(1-a)*(1+sin(linspace(-1,1,nx)*pi/2)).'/2;

[vf,pf] = vortFUnsteady(xf,solStruct);
[va,pa] = vortAUnsteady(xa,solStruct);

liftMat(j) = integral(@(xVar) vortFUnsteady(xVar.',solStruct).',-1,a) ...
            +integral(@(xVar) vortAUnsteady(xVar.',solStruct).', a,1);

xFMat(j,:) = xf;
xAMat(j,:) = xa;
pAMat(j,:)=pa;
pFMat(j,:)=pf;
vAMat(j,:)=va;
vFMat(j,:)=vf;

end


%% Plots
cols = hot(ceil(1.5*nJunct));
% figure(1)
% clf
% for j=1:nJunct
% plot(xAMat(j,:),abs(vAMat(j,:)),'LineWidth',2,'Color',cols(j,:));
% hold on
% ax1 = gca;
% ax1.ColorOrderIndex = j;
% plot(xFMat(j,:),abs(vFMat(j,:)),'LineWidth',2,'Color',cols(j,:));
% end
% hold off
% axis([-1,1,0,5])
% 
% ax2 = axes('Position',[.6 .6 .25 .25]);
% for j = 1:nJunct
%     scatter(ax2,aVec(j),abs(liftMat(j)), 'MarkerFaceColor',cols(j,:),'MarkerEdgeColor','k')
%      hold(ax2,'on');
% end
% 
% xlim([-1,1])
% xlabel(ax2,'$a$','Interpreter','latex')
% ylabel(ax2,'$\mathcal{L}$','Interpreter','latex')
% hold off
% 
% cleanfigure;
% matlab2tikz([imageFolder,'insetTest','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','axis on top');

figure(2)
clf
for j = 1:nJunct
plot(xAMat(j,:),abs(pAMat(j,:)),'LineWidth',2,'Color',cols(j,:));
hold on
ax1 = gca;
ax1.ColorOrderIndex = j;
plot(xFMat(j,:),abs(pFMat(j,:)),'LineWidth',2,'Color',cols(j,:));
end
hold off
axis([-1,1,.0,5])
xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$\Delta p$','Interpreter','latex')


ax2 = axes('Position',[.6 .6 .25 .25]);
for j = 1:nJunct
    plot(ax2,xAMat(j,:),psiFunCell{j}(xAMat(j,:)), '-', 'LineWidth', 2,'Color',cols(j,:))
     hold(ax2,'on');
    plot(ax2,xFMat(j,:),psiFunCell{j}(xFMat(j,:)), '-', 'LineWidth', 2,'Color',cols(j,:))
end
xlim([-1,1])
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$\psi$','Interpreter','latex')
hold off

cleanfigure;
% matlab2tikz([imageFolder,'insetTest',num2str(m),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','axis on top');

end

