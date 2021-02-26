function solStruct = calculateUnsteadyCoefficients(struct)

z = struct.z;
dzdx = struct.dzdx;

k = struct.k;
N = struct.N;

rhoe = struct.rhoe;
Phifun = struct.Phifun;
psiFun = @(x) eps + 4./(Phifun(x) + 2i*k.*rhoe(x));
psiFunQS =  @(x) eps + 4./(Phifun(x) + 0i*k.*rhoe(x));

wa = @(x) 1i*k*z(x)+dzdx(x);

fN = @(x)  2*wa(x);    
fC=  @(x)  1i*k*exp(1i*k*(1-x))/pi.*expint(-1i*k*(x-1));

if k<inf
leSing = 1/pi*acot(psiFun(-1)); % Leading edge singularity
teZero = 1/pi*acot(psiFun( 1)); % Trailing edge zero
else
   leSing = 0.5; teZero = 0.5; 
end
leSingQS = 1/pi*acot(psiFunQS(-1)); % Leading edge singularity
teZeroQS = 1/pi*acot(psiFunQS( 1)); % Trailing edge zero

small = 1e-10;

if k<inf

xCol = real(myJacobiNodes(N+1,teZero,-leSing)); % Generate collocation points

PGen  = myJacobiP(N+1,N-1,teZero,1-leSing,xCol);  
IGen  = myJacobiI(N-1,teZero,1-leSing,xCol);  
QGen  = myJacobiQ2(N-1,teZero,1-leSing,xCol);

PGenQS  = myJacobiP(N+1,N-1,teZeroQS,1-leSingQS,xCol);  
QGenQS  = myJacobiQ2(N-1,teZeroQS,1-leSingQS,xCol);

ISing = myJacobiI(1,teZero,-leSing,xCol);
QSing = myJacobiQ2(1,teZero,-leSing,xCol);

QSingQS = myJacobiQ2(1,teZeroQS,-leSingQS,xCol);

INC = myJacobiI(1,teZero-1,1-leSing,xCol);
QNC = myJacobiQ2(1,teZero-1,1-leSing,xCol);

IWake = myJacobiI(1,small,1-leSing,xCol);
QWake =  myJacobiQ2(1,small,1-leSing,xCol);

else
    
xCol = myJacobiNodes(N+1,0.5,-0.5); % Generate collocation points

end

%% Full problem

if k<inf

% Set up the collocation matrix
LwakeFull = -psiFun(xCol).*(weight(xCol,small,1-leSing) + 1i*k.*IWake(:,1)) ...
            + weight(xCol,small,1-leSing)/pi.*QWake(:,1);
LgenFull =  -psiFun(xCol).*(weight(xCol,teZero,1-leSing).*PGen + 1i*k.*IGen) ...
            + weight(xCol,teZero,1-leSing)/pi.*QGen;
LlesingFull = -psiFun(xCol).*(weight(xCol,teZero, -leSing) + 1i*k.*ISing(:,1)) ...
            + weight(xCol,teZero, -leSing)/pi.*QSing(:,1);

% Define the components of the wake coefficient        
Gamma0 = 2^(1+teZero-leSing)*double(beta(sym(1-leSing),sym(1+teZero)))/(1+2i*k*double(beta(sym(2-leSing),1)));
Gamma1 = 2^(2+teZero-leSing)*double(beta(sym(2-leSing),sym(1+teZero)))/(1+2i*k*double(beta(sym(2-leSing),1)));    
c2 = -1i*k*2^(leSing-1);    

Ffull = [LlesingFull     + Gamma0*(c2 * LwakeFull - fC(xCol)),...
         LgenFull(:,1)   + Gamma1*(c2 * LwakeFull - fC(xCol)),...
         LgenFull(:,2:end)];

% Solve to find the coefficients
fullCoefs  = Ffull\fN(xCol);
norm(Ffull*fullCoefs - fN(xCol),'inf');

% Calculate the coefficient of the wake
Gamma =  fullCoefs(1)*Gamma0 + fullCoefs(2)*Gamma1;

else
   Gamma = nan;
   fullCoefs = nan(size(xCol));
end

%% Non-circulatory problem

if k < inf
% Set up the collocation matrix
LteSingNC = -psiFun(xCol).*(weight(xCol,teZero-1,1-leSing) + 1i*k.*INC(:,1)) ...
            + weight(xCol,teZero-1,1-leSing)/pi.*QNC(:,1);
LgenNC =  LgenFull;
LlesingNC = LlesingFull;

% Define the components of the trailing-edge singularity coefficient        
denom = -2^(1+teZero-leSing)*double(beta(sym(teZero),sym(2-leSing)));
theta0 = 2^(1+teZero-leSing)*double(beta(sym(1-leSing),sym(1+teZero)))/denom;
theta1 = 2^(2+teZero-leSing)*double(beta(sym(2-leSing),sym(1+teZero)))/denom;

FNC = [LlesingNC     + theta0 * LteSingNC,...
       LgenNC(:,1)   + theta1 * LteSingNC,...
       LgenNC(:,2:end)];
   
% Solve to find the coefficients
ncCoefs  = FNC\fN(xCol);

elseif k == inf
    
IGen  = myJacobiI(N-1,0.5,0.5,xCol);  
QGen  = myJacobiQ2(N-1,0.5,0.5,xCol);

ISing = myJacobiI(1,0.5,-0.5,xCol);
QSing = myJacobiQ2(1,0.5,-0.5,xCol);

INC = myJacobiI(1,-0.5,0.5,xCol);
QNC = myJacobiQ2(1,-0.5,0.5,xCol);
    
% Set up the collocation matrix
if all(abs(Phifun(xCol))>1e-15)
LteSingNC = -2./rhoe(xCol).*INC(:,1) ...
            + weight(xCol,-0.5,0.5)/pi.*QNC(:,1);
LgenNC =  -2./rhoe(xCol).*IGen ...
            + weight(xCol,0.5,0.5)/pi.*QGen;
LlesingNC = -2./rhoe(xCol).*ISing(:,1) ...
            + weight(xCol,0.5,-0.5)/pi.*QSing(:,1);
else
LteSingNC =  weight(xCol,-0.5,0.5)/pi.*QNC(:,1);
LgenNC =   weight(xCol,0.5,0.5)/pi.*QGen;
LlesingNC = weight(xCol,0.5,-0.5)/pi.*QSing(:,1);    
end

% Define the components of the trailing-edge singularity coefficient        
denom = -2*double(beta(sym(0.5),sym(1.5)));
theta0 = (2*double(beta(sym(0.5),sym(1.5))))/denom;
theta1 = (2^2*double(beta(sym(1.5),sym(1.5))))/denom;    

FNC = [LlesingNC     + theta0 * LteSingNC,...
       LgenNC(:,1)   + theta1 * LteSingNC,...
       LgenNC(:,2:end)];
   
% Solve to find the coefficients
ncCoefs  = FNC\(2i*z(xCol));    
    
end

% Calculate the coeficient of the trailing-edge singularity
Theta = theta0 * ncCoefs(1) + theta1 * ncCoefs(2);

%% Quasi-static problem

if k<inf

% Set up the collocation matrix
LgenQS =  -psiFunQS(xCol).*weight(xCol,teZeroQS,1-leSingQS).*PGenQS ...
            + weight(xCol,teZeroQS,1-leSingQS)/pi.*QGenQS;
LsingQS = -psiFunQS(xCol).*weight(xCol,teZeroQS, -leSingQS) ...
            + weight(xCol,teZeroQS, -leSingQS)/pi.*QSingQS(:,1);

FQS = [LsingQS,LgenQS];

% Solve to find the coefficients
qsCoefs  = FQS\fN(xCol);
norm(FQS*qsCoefs - fN(xCol),'inf');
else
    qsCoefs = nan(size(xCol));
end
%% Put all the data into a structure
solStruct = struct;
solStruct.teZero = teZero;
solStruct.leSing = leSing;
solStruct.teZeroQS = teZeroQS;
solStruct.leSingQS = leSingQS;
solStruct.xCol = xCol;
solStruct.fullCoefs  = fullCoefs;
solStruct.qsCoefs  = qsCoefs;
solStruct.ncCoefs  = ncCoefs;
solStruct.Gamma = Gamma;
solStruct.Theta = Theta;
solStruct.psiFun = psiFun;


end