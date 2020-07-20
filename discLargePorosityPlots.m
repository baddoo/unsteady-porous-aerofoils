% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';

%% Calculate the p-coefficients
clf

struct.nf = 20; struct.na = 20;
nj = 5;
psiVec = logspace(-2,2,nj);
aVec = linspace(-1,1,nj+2);aVec(1) = []; aVec(end) = [];
nx = 100;
xAMat = zeros(nj,nx);xFMat = zeros(nj,nx);
pAMat = zeros(nj,nx);pFMat = zeros(nj,nx);
vAMat = zeros(nj,nx);vFMat = zeros(nj,nx);
psiFunCell = cell(1,nj);
a = 0;
for m = 1:2
    if m ==1; beta0 = 1; beta1 = beta0/2; elseif m==2; beta0 = 1; beta1 = 0; end
% m ==1 is pure pitching, m==2 is pure pitching.
z = @(x) beta0/2 + beta1*x; struct.z = z;
dzdx = @(x) beta1 + 0*x; struct.dzdx = dzdx;

for j = 1:nj

struct.junction = a;
if m ==1
    struct.k = .5;
elseif m==2
    struct.k = .1;
end
psiFun = @(xVar) heaviside(xVar-a).*psiVec(j) + eps;
psiFunCell{j} = psiFun;

rhoe = @(x) 1.2 + 0*x;
struct.Phifun = @(x) 1./psiFun(x);
struct.rhoe = rhoe;


solStruct = calculateUnsteadyCoefficientsDiscont(struct);   

xf = -1+(1+a)*(1+sin(linspace(-1,1,nx)*pi/2)).'/2;
xa = a+(1-a)*(1+sin(linspace(-1,1,nx)*pi/2)).'/2;

[vf,pf] = vortFUnsteady(xf,solStruct);
[va,pa] = vortAUnsteady(xa,solStruct);

xFMat(j,:) = xf;
xAMat(j,:) = xa;
pAMat(j,:)=pa;
pFMat(j,:)=pf;
vAMat(j,:)=va;
vFMat(j,:)=vf;

end

%% Rigid cases

rigStruct = struct;
rigStruct.type = 'full'; rigStruct.N = 40;
rigStruct.Phifun = @(xVar) 1./eps;
rigStruct2 = rigStruct;
f = @(xVar) .5*(a*(xVar+1)+xVar-1);
finv = @(xVar) (2*xVar+1-a)/(a + 1);
rigStruct2. k = (struct.k)/2*(a+1); rigStruct2.z = @(xVar) 2/(a+1)*(z(f(xVar)));
rigStruct2.dzdx = @(xVar) dzdx(f(xVar));
rigSolStruct = calculateUnsteadyCoefficients(rigStruct); 
rigSolStruct2= calculateUnsteadyCoefficients(rigStruct2); 

xRig = sin(linspace(-1,1)*pi/2); xRig(1) = []; xRig(end) = [];
rigSol = presFun(xRig',rigSolStruct)';
rigSol2= presFun(xRig',rigSolStruct2)';

%% Plots
cols = hot(ceil(1.5*nj));

figure(m)
clf
hold on
plot(xRig,abs(rigSol),'k','LineWidth',2)
plot(f(xRig),abs(rigSol2),'k','LineWidth',3)

for j = 1:nj
%    yyaxis left
plot(xAMat(j,:),abs(pAMat(j,:)),'-','LineWidth',1,'Color',cols(j,:));
ax1 = gca;
ax1.ColorOrderIndex = j;
plot(xFMat(j,:),abs(pFMat(j,:)),'-','LineWidth',1,'Color',cols(j,:));
%yyaxis right
%plot(xFMat(j,:),psiFunCell{j}(xFMat(j,:)).*abs(pFMat(j,:)), '--', 'LineWidth', 2,'Color',cols(j,:))
%plot(xAMat(j,:),psiFunCell{j}(xAMat(j,:)).*abs(pAMat(j,:)), '--', 'LineWidth', 2,'Color',cols(j,:))
end
hold off
xlim([-1,1]);
%yyaxis left
if m==1; ylim([0,5]); elseif m ==2; ylim([0,.5]); end
%yyaxis right    
%if m==1; ylim([0,.4]); elseif m ==2; ylim([0,40]); end

xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$|\Delta p|$','Interpreter','latex')
grid on

ax2 = axes('Position',[.6 .6 .25 .25]);
for j = 1:nj
    semilogy(ax2,xAMat(j,2:end),psiFunCell{j}(xAMat(j,2:end)), '-', 'LineWidth', 2,'Color',cols(j,:))
     hold(ax2,'on');
    %semilogy(ax2,xFMat(j,:),psiFunCell{j}(xFMat(j,:)), '-', 'LineWidth', 2,'Color',cols(j,:))
end
xlim([-1,1])
ylim([.005,120])
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$1/\Phi$','Interpreter','latex')
hold off
grid minor

 cleanfigure;
 matlab2tikz([imageFolder,'largePor',num2str(m),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
                  'extraaxisoptions','axis on top=false');

end

