% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi/images/';

% Set angle of attack and parabolic camber (only use alp=0 for now)
beta0 = 1; beta1 = 0;

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) beta1 + 0*xVar; struct.dzdx = dzdx;

kVec = [.1,1,3];
%% Calculate the p-coefficients
clf

nJunctI = 5;
aVecI = linspace(-1,1,nJunctI+2);aVecI(1) = []; aVecI(end) = [];
aVecL = linspace(-.99,.99,100);
aVec = unique(sort([aVecI,aVecL]));
locs = zeros(1,nJunctI); for j = 1:nJunctI, locs(j) = find(aVecI(j)==aVec); end
nJunctF = numel(aVec);
nx = 100;
xAMat = zeros(nJunctI,nx);xFMat = zeros(nJunctI,nx);
pAMat = zeros(nJunctI,nx);pFMat = zeros(nJunctI,nx);
vAMat = zeros(nJunctI,nx);vFMat = zeros(nJunctI,nx);
liftMat = zeros(3,nJunctF);
psiFunCell = cell(3,nJunctI);

for m = [1]%,2,3]
    nk = 1;
for j = 1:nJunctF
struct.k = kVec(m);

a = aVec(j);
struct.nf =ceil((1+a)/2*20)+1; struct.na =ceil((1-a)/2*20)+1;

struct.junction = a;
if m ==1
struct.psiFun = @(xVar) .1i+.15*heaviside(xVar-a)./(1-a);
elseif m==2
struct.psiFun = @(xVar) eps+.1*heaviside(xVar-a).*(1/(1-a)+1.2*(xVar-a)/(1-a).^2/2);
elseif m==3
struct.psiFun = @(xVar) eps+.1*((1+xVar)+heaviside(xVar-a)/(1-a));
end
solStruct = calculateUnsteadyCoefficientsDiscont(struct);   

liftMat(m,j) = integral(@(xVar) vortFUnsteady(xVar.',solStruct).',-1,a) ...
                +integral(@(xVar) vortAUnsteady(xVar.',solStruct).', a,1);

if ismember(j,locs)
    psiFunCell{m,nk} = struct.psiFun;

    xf = [-1+(1+a)*(1+sin(linspace(-1,1,nx-1)*pi/2))/2,a-eps].';
    xa = [a+eps,a+(1-a)*(1+sin(linspace(-1,1,nx-1)*pi/2))/2].';

    pf = presFUnsteady(xf,solStruct);
    pa = presAUnsteady(xa,solStruct);

    xFMat(nk,:) = xf;
    xAMat(nk,:) = xa;
    pAMat(nk,:)=pa;
    pFMat(nk,:)=pf;
        nk = 1 + nk;

end

end


%% Plots
if m ==1
cols = hot(ceil(1.5*nJunctI));
elseif m==2
cols = winter(ceil(nJunctI));
elseif m==3
cols = pink(ceil(1.5*nJunctI));
end
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
for j = 1:nJunctI
plot(xAMat(j,:),abs(pAMat(j,:)),'LineWidth',2,'Color',cols(j,:));
hold on
ax1 = gca;
ax1.ColorOrderIndex = j;
plot(xFMat(j,:),abs(pFMat(j,:)),'LineWidth',2,'Color',cols(j,:));
end
hold off
xlim([-1,1])
if m==1; ylim([0,.75]);elseif m==2; ylim([0,5]); elseif m==3; ylim([0,35]); end

xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$\Delta p$','Interpreter','latex')

ax2 = axes('Position',[.4 .6 .45 .25]);
for j = 1:nJunctI
    plot(ax2,xAMat(j,:),psiFunCell{m,j}(xAMat(j,:)), '-', 'LineWidth', 2,'Color',cols(j,:))
     hold(ax2,'on');
    plot(ax2,xFMat(j,:),psiFunCell{m,j}(xFMat(j,:)), '-', 'LineWidth', 2,'Color',cols(j,:))
end
xlim([-1,1])
ylim([0,.6])
%if m==1; ylim([0,.33]); elseif m==2; ylim([0,.2]); end
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$\psi$','Interpreter','latex')
hold off

cleanfigure;
%matlab2tikz([imageFolder,'insetTest',num2str(m),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                 'extraaxisoptions','axis on top');

end

%% Lift plot
figure(3)
normFun = @(n,mat) abs((mat(n,:)-mat(n,(end+1)/2))./mat(n,(end+1)/2));
normFun = @(n,mat) abs(mat(n,:));
plot(aVec,normFun(1,liftMat),'r','LineWidth',2)
hold on
plot(aVec,normFun(2,liftMat),'b','LineWidth',2)
colPink = pink;
plot(aVec,normFun(3,liftMat),'LineWidth',2,'Color',colPink(20,:))
hold off

xlim([-1,1])
%ylim([0,.33])
xlabel('$c$','Interpreter','latex')
ylabel('$\Delta \mathcal{L}$','Interpreter','latex')

cleanfigure;
%matlab2tikz([imageFolder,'constPorLift.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                 'extraaxisoptions','axis on top');
