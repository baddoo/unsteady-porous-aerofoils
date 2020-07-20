
beta0 = 1; beta1 = 0;

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) (beta1 + 0*xVar); struct.dzdx = dzdx;

nVar = 10;

struct.N = 40;

k = 1;
struct.k = k;
a = .5;
psiFunCell = @(xVar) eps + 0*heaviside(real(xVar-a)).*((xVar-a).^2);
struct.psiFun = psiFunCell;

qsStruct = struct;
qsStruct.type = 'nc';
ncSolStruct = calculateUnsteadyCoefficientsNC(qsStruct); 

tic
qsLift= quadgk(@(xVar) vortFun(xVar.',ncSolStruct).',-1,1,'AbsTol',1e-15)
%qsLift= quadgk(@(xVar) vortFun(xVar.',solStructQS).',-1,1)
%qsLift= integral(@(xVar) vortFun(xVar.',solStructQS).',-1,1)
toc
tic
circulation(ncSolStruct)
toc
return
ncSolStruct= calculateUnsteadyCoefficientsNC(struct);   
ncLift = quadgk(@(xVar) presFunNC(xVar.',ncSolStruct).',-1,1,'AbsTol',1e-12);

struct.type = 'full';
fullSolStruct = calculateUnsteadyCoefficients(struct); 
fullLift = quadgk(@(xVar) presFun(xVar.',fullSolStruct).',-1,1,'AbsTol',1e-12);
