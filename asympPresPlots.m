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
beta0 = 1; beta1 = 2; beta2 = 3;
z = @(xVar) beta0/2 + beta1*xVar + beta2*xVar.^2; struct.z = z;
dzdx = @(xVar) beta1 + 0*xVar + 2*beta2*xVar; struct.dzdx = dzdx;
nk = 1;
fullCirc = zeros(1,nk);
ncCirc = zeros(1,nk);
qsCirc = zeros(1,nk);
fullLift = zeros(1,nk);
ncLift = zeros(1,nk);
qsLift = zeros(1,nk);
kVec = logspace(0,2,nk);
kVec = 50;

for j = 1:nk
    k = kVec(j);
    struct.k = k;
    struct.N = round(10+20*sqrt(k)); 
rhoe = @(x) 1 + 0*x;
PSIfun = @(x) 1-x;
struct.PSIfun = PSIfun;
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
%xp = 1-2*logspace(-5,0,30)';
A = 1i*infSolStruct.Theta*sqrt(2)*k.^1.5;
omega = 1i*A/k*sqrt(pi).*exp(-1i*(k+pi/4));
[fullVInf,ncVInf] = vortFun(xp,infSolStruct);
xe = (xp-1)*k;
outerV = ncVInf.*k;
innerV =  -1i*k*omega*exp(1i*(k-xe)).*(1 - double(erf(sym(exp(-1i*pi/4)*1i*sqrt(-xe)))));
S = A./sqrt(xe);
compV =  outerV + 0*innerV - 0*S;

[~,ncPInf] = presFun(xp,infSolStruct);

[~,ncInt,~] = vortInt(xp,infSolStruct);

compP = -2*double(compV + 1i*k*(k*ncInt + 0*omega*exp(1i*k).*...
            (exp(-1i*xe).*(1 - erf((sym(exp(-1i*pi/4)*1i*sqrt(-xe)))))...
            -exp( 2i*k ).*(1 - erf((sym(exp(-1i*pi/4)*1i*sqrt(--2*k))))))));
        
figure(2)
[fullP,ncP,qsP] = presFun(xp,solStruct);
[fullV,ncV,qsV] = vortFun(xp,solStruct);

figure(2)
%plot(xp,real(fullV)/k,'r');
%semilogy(xp,abs(fullV).*(1-xp.^2),'LineWidth',3);
%hold on
plot(xp,abs(fullP-0*k^2*ncPInf)./k^2,'r--','LineWidth',3);
hold on
%plot(xp,real(ncV)/k,'b');
%plot(xp,imag(ncV)/k,'b--');
plot(xp,abs(compP-0*k^2*ncPInf)./k^2,'g')
%plot(xp,abs(ncP)/k^2,'b')
%plot(xp,imag(outerSol./k),'b')
%plot(xp,imag(innerSol./k),'s')
%semilogy(xp,abs(ncV),'LineWidth',3);
hold off
xlim([-1,1])
%ylim([-1,10])

%%
return
figure(3)
plot(kVec,abs(fullLift-ncLift)./kVec)
hold on
%plot(kVec,real(ncLift)./kVec.^2
plot(kVec,abs(qsLift)./kVec)
hold off