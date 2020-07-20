% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion

% Set angle of attack and parabolic camber (only use alp=0 for now)
beta0 = 0; beta1 = 1;

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) beta1 + 0*xVar; struct.dzdx = dzdx;

%% Calculate the p-coefficients
clf

nVar = 100;

solCell = cell(1,nVar);

struct.N = 50;
nPsi = 1;
profile on
tic
freqVec = linspace(0,5,nVar+1); freqVec(1) = [];
Uvec = 1./freqVec;
numLES = zeros(nVar,nPsi);
numPres = zeros(nVar,nPsi);
numPow = zeros(nVar,nPsi);

psi = eps;%linspace(0,2,nPsi);
for m = 1:nPsi
for l = 1:nVar

k=2*pi*freqVec(l);
struct.k = k;
struct.psiFun = @(xVar) eps+ psi(m)*(1+xVar)/2;
solStruct = calculateUnsteadyCoefficients(struct);   
solCell{l} = solStruct;
%fullCell{l} = @(zVar) presFun(zVar.',solStruct).';
numLES(l,m) = pi^3./k^2*(real(solStruct.coefs(1)).^2 + imag(solStruct.coefs(1)).^2);
presInt = @(xVar) real(presFun(xVar.',solStruct)).'.*real(dzdx(xVar)) ...
                 + imag(presFun(xVar.',solStruct)).'.*imag(dzdx(xVar));
powInt = @(xVar) real(presFun(xVar.',solStruct)).'.*real(1i*k*z(xVar)) ...
                 + imag(presFun(xVar.',solStruct)).'.*imag(1i*k*z(xVar));
numPres(l,m) = 2*pi^2./k^2*integral(presInt,-1,1);
numPow(l,m) = 2*pi^2./k^2*integral(powInt,-1,1);

end
end
toc
profile off

sig = 2*pi./Uvec;
C = @(sigVar) besselk(1,1i*sigVar)./(besselk(0,1i*sigVar) + besselk(1,1i*sigVar));

a = zeros(3,nVar);
a(1,:) = -2i*pi*Uvec.*C(sig)*beta0 + 2i*pi*Uvec.*(1-C(sig))*beta1 - 2*C(sig).*Uvec.^2*beta1;
a(2,:) = 2*pi^2*beta0 - 4i*pi*Uvec*beta1;
a(3,:) = pi^2*beta1;

mooreLES = pi./(4*Uvec.^2).*(real(a(1,:)).^2 + imag(a(1,:)).^2);
moorePres = pi/2*(real(a(1,:)).*real(beta1) + imag(a(1,:)).*imag(beta1))...
           +pi^3*(real(beta0)*real(beta1)+imag(beta0)*imag(beta1));
mooreThrust = mooreLES + moorePres;
numThrust = numLES + numPres;

%%
cols = hot(ceil(1.5*(nPsi)));
semilogy(freqVec,numThrust(:,1),'k','LineWidth',2)
%semilogy(freqVec,-numThrust(:,1)./numPow(:,1).*Uvec.','k','LineWidth',2)
hold on
for l = 1:nPsi
%semilogy(freqVec,-numThrust(:,l)./numPow(:,l).*Uvec.','Color',cols(l,:),'LineWidth',2)
end
semilogy(freqVec,abs(mooreThrust),'r--','LineWidth',2)
hold off

%%
return
figure(2) 
error = zeros(nVar-1,1);
for l0 = 1:(nVar-1)
     intFun = @(xVar) abs(fullCell{l0+1}(xVar) - fullCell{l0}(xVar)).^2;
     norFun = @(xVar) abs(fullCell{l0+1}(xVar)).^2;
     error(l0)= sqrt(quadgk(intFun,-1,1)./quadgk(norFun,-1,1));
    %error(l0) = sqrt(max(intFun(xPlot))./max(norFun(xPlot)));
end

length = 1:nVar-1;
loglog(length,abs(error));
