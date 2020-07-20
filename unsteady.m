% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion

% Set angle of attack and parabolic camber (only use alp=0 for now)
aoa = 1; cam = 0;
k=2;

% Set porosities for evaluation

% Defines geometry and derivative
zGeo = @(xVar) -cam/2*(xVar.^2-1);
z = @(xVar) -aoa*xVar + cam*(1-xVar.^2)/2;
dzdx = @(xVar) -aoa - cam*xVar;
wa = @(xVar) 1i*k*z(xVar)+dzdx(xVar);

%% Calculate the p-coefficients
nVar = 1;

figure(6)
clf

psiV = 0*logspace(-2,0,nVar);
psiFunCell = cell(1,nVar);
%kVec = linspace(0,5,nVar);
for l = 1:nVar
    psiFunCell{l} = @(xVar) 2+psiV(l) + 0*xVar;
    %psiFunCell{l} = @(xVar) psiV(l)*(1+xVar)/2;
    %psiFunCell{l} = @(xVar) psiV(l).*(1+xVar).^2/4;
end

cFunCell = cell(1,nVar);
ncFunCell = cell(1,nVar);
fullCell = cell(1,nVar);

wakeStrength = zeros(1,nVar);
N = 25;

profile on
tic
for l = 1:nVar
fN = @(xVar) -2*wa(xVar);    
fNC= @(xVar)  -1i*k*exp(-1i*k*xVar)/pi.*expint(-1i*k*(xVar-1));
psi = psiFunCell{l};

leSing = 1/pi*acot(psi(-1)); % Leading edge singularity
teZero = 1/pi*acot(psi( 1)); % Trailing edge zero

const = -1i*k*2^leSing*exp(-1i*k);

xCol = myJacobiNodes(N+1,teZero,1-leSing); % Generate collocation points

jP  = myJacobiP(N+1,N-1,teZero,1-leSing,xCol);  
jI  = myJacobiI(N-1,teZero,1-leSing,xCol);  
jIC  = myJacobiI(1,teZero-1,-leSing,xCol);  
jQ  = myJacobiQ2(N-1,teZero,1-leSing,xCol);
jQC  = myJacobiQ2(1,teZero-1,-leSing,xCol);
jPVar = myJacobiP(N+1,1,teZero,-leSing,xCol);
jIVar = myJacobiI(1,teZero,-leSing,xCol);
jIop = myJacobiI(1,0,-leSing,xCol);
jQVar = myJacobiQ2(1,teZero,-leSing,xCol);
jQop =  myJacobiQ2(1,1e-5,-leSing,xCol); % address the fact that we can't use zero parameter

% Make this into a nice function
op =    psi(xCol).*(weight(xCol,1e-5,-leSing)        +1i*k.*jIop(:,1)) ...
        - weight(xCol,1e-5,-leSing)/pi.*jQop(:,1);
fC =    psi(xCol).*(weight(xCol,teZero-1,-leSing)        +1i*k.*jIC(:,1)) ...
        - weight(xCol,teZero-1,-leSing)/pi.*jQC(:,1);
f =     psi(xCol).*(weight(xCol,teZero,1-leSing).*jP      +1i*k.*jI) ...
        - weight(xCol,teZero,1-leSing)/pi.*jQ;
f0 = psi(xCol).*(weight(xCol,teZero, -leSing).*jPVar(:,1) +1i*k.*jIVar(:,1)) ...
        - weight(xCol,teZero, -leSing)/pi.*jQVar(:,1);
    
F = [f0    - 2*(teZero./(1+teZero-leSing)).*fC,...
    f(:,1) - 4*teZero*(1-leSing)./(1+teZero-leSing)./(2+teZero-leSing).*fC,...
    f(:,2:end)];

lam1 = const*2^(1-teZero)*beta(1-leSing,1)/beta(1-leSing,teZero) + exp(-1i*k)/(2^(teZero-leSing)*beta(1-leSing,teZero));

ncCoefs  = F\ fN(xCol);
cCoefs   = F\(fNC(xCol) - const.*op - lam1*fC);

lamNC= -2^(-leSing)*(ncCoefs(1)*2*teZero/(1+teZero-leSing)...
           +ncCoefs(2)*4*teZero*(1-leSing)/(1+teZero-leSing)/(2+teZero-leSing));
lamC=  -2^(-leSing)*(cCoefs(1)*2*teZero/(1+teZero-leSing)...
           + cCoefs(2)*4*teZero*(1-leSing)/(1+teZero-leSing)/(2+teZero-leSing)...
           +(2*const*beta(1-leSing,1)-2^leSing*exp(-1i*k))/(2^teZero*beta(1-leSing,teZero)));
       
wakeStrength(l) = -lamNC/lamC;

ncFunL = @(zVar) ncFun(zVar.',ncCoefs, teZero,leSing).';
cFunL =  @(zVar)  cFun(zVar.', cCoefs, teZero,leSing,const,k).';
ncFunCell{l}  = ncFunL;
cFunCell{l}  = cFunL;
fullCell{l} = @(zVar) ncFunL(zVar) + wakeStrength(l)*cFunL(zVar);
   %abs(quadgk(ncFunCell{l},-1,1))
   %(integral(cFunCell{l},-1,1,'RelTol',1e-12,'AbsTol',1e-12))

end
toc
profile off

wakeInt = @(xVar) sqrt((1+xVar)./(1-xVar)).*wa(xVar);
rigidWakeStrength = 4*quadgk(wakeInt,-1,1)/(1i*pi*k*(besselh(1,2,k)+1i*besselh(0,2,k)));
abs(wakeStrength(l))./abs(rigidWakeStrength)

%% Plots
nP = 500;
xPlot = sin(pi/2*linspace(-1,1,nP+2)); xPlot(1) = []; xPlot(end) = [];
cols = (hot(ceil(1.5*nVar)));

clf; cla

ax1 = gca;
for l0 = 1:nVar
    %ncvals = ncFunCell{l0}(xPlot);
    %cvals = cFunCell{l0}(xPlot);
    %semilogy(ax1,xPlot, abs(ncvals + wakeStrength(l0)*cvals),'-' ,'LineWidth',4,'Color',cols(l0,:))    
    plot(ax1,xPlot, real(fullCell{l0}(xPlot)),'-' ,'LineWidth',4,'Color',cols(l0,:))    

    hold(ax1,'on');
end
%set(ax1,'xlim',[.99,1])
set(ax1,'ylim',[-1,2])
xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$\gamma$','Interpreter','latex')

ax2 = axes('Position',[.6 .6 .25 .25]);
for l0 = 1:nVar
    plot(ax2,xPlot,psiFunCell{l0}(xPlot) , '-', 'LineWidth', 2,'Color',cols(l0,:))
     hold(ax2,'on');
end
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$\psi$','Interpreter','latex')
hold off
%%
return

figure(2) 
error = zeros(nVar-1,1);
for l0 = 1:(nVar-1)
     intFun = @(xVar) abs(fullCell{l0+1}(xVar) - fullCell{l0}(xVar)).^1;
     norFun = @(xVar) abs(fullCell{l0+1}(xVar)).^1;
     error(l0)= sqrt(quadgk(intFun,-1,1)./quadgk(norFun,-1,1));
    %error(l0) = sqrt(max(intFun(xPlot))./max(norFun(xPlot)));
end

length = 1:nVar-1;
semilogy(length,abs(error));