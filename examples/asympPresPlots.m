% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi/images/';

%% Calculate the p-coefficients
clf

nPsi = 1;
psi = linspace(0,2,nPsi);
N = 30;
struct.N = N; 
struct.type = 'full';
nj = 5;
PSIVec = linspace(0,5,nj);
xMat = zeros(nj,N+1);
pMat = zeros(nj,N+1);
vMat = zeros(nj,N+1);
psiFunCell = cell(1,nj);
a = 0;

profile on

    
% Set angle of attack and parabolic camber (only use alp=0 for now)
beta0 = 1; beta1 = 1; beta2 = 2;
z = @(xVar) beta0/2 + beta1*xVar + beta2*xVar.^2; struct.z = z;
dzdx = @(xVar) beta1 + 2*beta2*xVar; struct.dzdx = dzdx;
nk = 1;
fullCirc = zeros(1,nk);
ncCirc = zeros(1,nk);
qsCirc = zeros(1,nk);
fullLift = zeros(1,nk);
ncLift = zeros(1,nk);
qsLift = zeros(1,nk);
kVec = logspace(0,2,nk);
kVec = 55;

for j = 1:nk
    k = kVec(j);
    struct.k = k;
    struct.N = round(10+20*sqrt(k)); 
rhoe = @(x) 1 + 0*x;
Phifun = @(x) 1+0*x;%./(1-x);
struct.Phifun = Phifun;
struct.rhoe = rhoe;
solStruct = calculateUnsteadyCoefficients(struct); 
[fCirc,nCirc,qCirc] = circulation(solStruct);
[fLift,nLift,qLift] = lift(solStruct);
fullCirc(j) = fCirc;
ncCirc(j) = nCirc;
qsCirc(j) = qCirc;
fullLift(j) = fLift;
ncLift(j) = nLift;
qsLift(j) = qLift;
j
end

% %%
% figure(1)
% subplot(1,2,1)
% loglog(kVec,abs(fullCirc),'b','LineWidth',3); 
% hold on; 
% loglog(kVec,(kVec/kVec(end)).^.5*abs(fullCirc(end)),'b--','LineWidth',3);
% %loglog(kVec,kVec.^.25,'b--','LineWidth',3);
% loglog(kVec,abs(qsCirc),'r','LineWidth',3);
% loglog(kVec,kVec/kVec(end)*abs(qsCirc(end)),'r--','LineWidth',3); 
% hold off
% title('circulation')
% 
% subplot(1,2,2)
% loglog(kVec,abs(fullLift),'b','LineWidth',3); 
% hold on; 
% loglog(kVec,kVec.^2,'b--','LineWidth',3); 
% loglog(kVec,abs(ncLift),'g','LineWidth',3); 
% loglog(kVec,kVec.^2,'g--','LineWidth',3); 
% loglog(kVec,abs(qsLift),'r','LineWidth',3); 
% loglog(kVec,kVec,'r--','LineWidth',3); 
% title('lift')
% 
% hold off

%%
struct.k = inf;
infSolStruct = calculateUnsteadyCoefficients(struct);

%%
xp = cos(flip(linspace(0,pi,50))).'; xp(1) = []; xp(end) = [];
[fullP,ncP,qsP] = presFun(xp,solStruct);
[fullV,ncV,qsV] = vortFun(xp,solStruct);

%xp = 1-2*logspace(-5,0,30)';
A = infSolStruct.Theta*sqrt(2);
[fullVInf,ncVInf] = vortFun(xp,infSolStruct);
xe = (xp-1)*k;
outerV = ncVInf;
circ = 1i*A*exp(-1i*pi/4)*sqrt(pi/k);
innerV = -1i*k*circ.*exp(1i*xe).*erfc(sym(exp(-1i*pi/4)*sqrt(-xe)));
S = A*sqrt(k./-xe);
%compV =  outerV + (innerV - S);
compV = outerV - 1i*k*circ*(-exp(1i*pi/4)./sqrt(pi*-xe)+ exp(1i*xe).*erfc(sym(exp(-1i*pi/4)*sqrt(-xe))));
[~,ncPInf] = presFun(xp,infSolStruct);
outerP = ncPInf;
innerAndOverlapP = -2i*(-1i*circ*-1i*exp(1i*xe).*erfc(sym(exp(-1i*pi/4)*sqrt(-xe))));
compP = outerP + innerAndOverlapP;

figure(1)
LW = 'LineWidth';
plot(xp,real(1i*fullV)/k,'k',LW,3);
hold on
plot(xp,real(1i*compV),'r',LW,2);
plot(xp,real(1i*outerV),'g--',LW,2);
plot(xp,real(1i*innerV),'m--',LW,2);
hold off

figure(2)
LW = 'LineWidth';
plot(xp,real(fullP)/k^2,'k',LW,3);
hold on
plot(xp,real(compP),'r',LW,2);
plot(xp,real(outerP),'g--',LW,2);
hold off
%plot(xp,real(1i*innerP),'m--',LW,2);


